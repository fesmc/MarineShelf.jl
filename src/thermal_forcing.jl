# ============================================================
#  thermal_forcing.jl
#  Freezing-point calculation, thermal-forcing assembly,
#  and basin-mean averaging for quad-nl / pico.
# ============================================================

# ------------------------------------------------------------
# Freezing point
# ------------------------------------------------------------

"""
    freezing_point(S, z_base, fp::FreezingPointParams{T}) → T

Pressure- and salinity-dependent freezing point [K] following
Favier et al. (2019), Eq. 3:

    T_fp = λ1·S + λ2 + λ3·z_base + T0

- `S`      — salinity [PSU]
- `z_base` — ice-shelf base elevation relative to sea level [m]
             (negative below sea level)
"""
@inline function freezing_point(
    S::T, z_base::T, fp::FreezingPointParams{T},
) where T<:AbstractFloat
    return fp.λ1 * S + fp.λ2 + fp.λ3 * z_base + T(T0)
end

"""
    freezing_point!(T_fp, S_shlf, z_base, fp)

In-place broadcast version over 2-D arrays.
"""
function freezing_point!(
    T_fp  ::AbstractMatrix{T},
    S_shlf::AbstractMatrix{T},
    z_base::AbstractMatrix{T},
    fp    ::FreezingPointParams{T},
) where T<:AbstractFloat
    @inbounds for j in axes(T_fp, 2), i in axes(T_fp, 1)
        T_fp[i,j] = freezing_point(S_shlf[i,j], z_base[i,j], fp)
    end
    return T_fp
end

# ------------------------------------------------------------
# Thermal forcing (absolute method)
# ------------------------------------------------------------

"""
    thermal_forcing_abs!(tf, T_shlf, T_fp, tf_corr, tf_corr_basin)

Thermal forcing for absolute-temperature methods (lin, quad, quad-nl):

    TF = T_shlf - T_fp + tf_corr + tf_corr_basin
"""
function thermal_forcing_abs!(
    tf           ::AbstractMatrix{T},
    T_shlf       ::AbstractMatrix{T},
    T_fp         ::AbstractMatrix{T},
    tf_corr      ::AbstractMatrix{T},
    tf_corr_basin::AbstractMatrix{T},
) where T<:AbstractFloat
    @. tf = T_shlf - T_fp + tf_corr + tf_corr_basin
    return tf
end

"""
    thermal_forcing_anom!(tf, dT_shlf, tf_corr, tf_corr_basin)

Thermal forcing for the anomaly method:

    TF = dT_shlf + tf_corr + tf_corr_basin
"""
function thermal_forcing_anom!(
    tf           ::AbstractMatrix{T},
    dT_shlf      ::AbstractMatrix{T},
    tf_corr      ::AbstractMatrix{T},
    tf_corr_basin::AbstractMatrix{T},
) where T<:AbstractFloat
    @. tf = dT_shlf + tf_corr + tf_corr_basin
    return tf
end

# ------------------------------------------------------------
# Basin-mean averaging
# ------------------------------------------------------------

"""
    basin_mean!(var_basin, var2D, f_grnd, basins, H_ice, mask_ocn)

Replace every point in each basin with the mean of `var2D` over
floating-ice points in that basin.

Averaging mask: `f_grnd < 1` AND (`H_ice > 0` OR point is
`FloatingLine` / `GroundingLine`).

If a basin has no qualifying points a warning is printed and the
basin mean is set to zero (this case is rare in practice).
"""
function basin_mean!(
    var_basin::AbstractMatrix{T},
    var2D    ::AbstractMatrix{T},
    f_grnd   ::AbstractMatrix{T},
    basins   ::AbstractMatrix{T},
    H_ice    ::AbstractMatrix{T},
    mask_ocn ::AbstractMatrix{OceanMask},
) where T<:AbstractFloat
    nx, ny = size(var2D)
    nb = Int(maximum(basins))

    for m in 1:nb
        s    = zero(T)
        npts = 0
        @inbounds for j in 1:ny, i in 1:nx
            if Int(basins[i,j]) == m &&
               f_grnd[i,j] < one(T) &&
               (H_ice[i,j] > zero(T) ||
                mask_ocn[i,j] == FloatingLine ||
                mask_ocn[i,j] == GroundingLine)
                s    += var2D[i,j]
                npts += 1
            end
        end

        mean_val = if npts > 0
            s / T(npts)
        else
            @warn "basin_mean!: no floating-ice points in basin $m; setting mean to 0"
            zero(T)
        end

        @inbounds for j in 1:ny, i in 1:nx
            Int(basins[i,j]) == m && (var_basin[i,j] = mean_val)
        end
    end
    return var_basin
end

"""
    basin_mean_icefree!(var_basin, var2D, f_grnd, basins, H_ice, mask_ocn)

Variant used by PICO: average over ice-free ocean points (`f_grnd == 0`
AND `H_ice == 0`) rather than floating-ice points.
"""
function basin_mean_icefree!(
    var_basin::AbstractMatrix{T},
    var2D    ::AbstractMatrix{T},
    f_grnd   ::AbstractMatrix{T},
    basins   ::AbstractMatrix{T},
    H_ice    ::AbstractMatrix{T},
    mask_ocn ::AbstractMatrix{OceanMask},
) where T<:AbstractFloat
    nx, ny = size(var2D)
    nb = Int(maximum(basins))

    for m in 1:nb
        s    = zero(T)
        npts = 0
        @inbounds for j in 1:ny, i in 1:nx
            if Int(basins[i,j]) == m &&
               f_grnd[i,j] == zero(T) &&
               H_ice[i,j]  == zero(T)
                s    += var2D[i,j]
                npts += 1
            end
        end

        mean_val = if npts > 0
            s / T(npts)
        else
            @warn "basin_mean_icefree!: no ice-free points in basin $m; setting mean to 0"
            zero(T)
        end

        @inbounds for j in 1:ny, i in 1:nx
            Int(basins[i,j]) == m && (var_basin[i,j] = mean_val)
        end
    end
    return var_basin
end

"""
    apply_basin_correction!(field, basins, corr_dict)

Add per-basin scalar corrections from `corr_dict :: Dict{Int,T}`
to `field` in-place.  Basins not present in the dict are unchanged.
"""
function apply_basin_correction!(
    field    ::AbstractMatrix{T},
    basins   ::AbstractMatrix{T},
    corr_dict::Dict{Int,T},
) where T<:AbstractFloat
    isempty(corr_dict) && return field
    @inbounds for j in axes(field, 2), i in axes(field, 1)
        b = Int(basins[i,j])
        if haskey(corr_dict, b)
            field[i,j] += corr_dict[b]
        end
    end
    return field
end
