# ============================================================
#  vertical_interp.jl
#  Vertical weighting and interpolation of ocean profiles
#  to a target depth (shelf-base or bed depth).
# ============================================================

# ------------------------------------------------------------
# Weight computation  (three methods, matching the Fortran)
# ------------------------------------------------------------

"""
    shelf_weights_mean!(wt, depth; depth_min, depth_max)

Uniform weighting of all depth layers within `[depth_min, depth_max]`.
`depth` is positive downward.  Returns normalised weights.
"""
function shelf_weights_mean!(
    wt       ::AbstractVector{T},
    depth    ::AbstractVector{T};
    depth_min::T,
    depth_max::T,
) where T<:AbstractFloat
    fill!(wt, zero(T))
    k0 = argmin(abs.(depth .- depth_min))
    k1 = argmin(abs.(depth .- depth_max))
    if k1 < k0
        error("shelf_weights_mean!: depth_min=$depth_min > depth_max=$depth_max " *
              "after nearest-layer matching (k0=$k0, k1=$k1).")
    end
    wt[k0:k1] .= one(T)
    return wt
end

"""
    shelf_weights_layer!(wt, depth, depth_shlf)

All weight on the single depth layer nearest to `depth_shlf`.
"""
function shelf_weights_layer!(
    wt         ::AbstractVector{T},
    depth      ::AbstractVector{T},
    depth_shlf ::T,
) where T<:AbstractFloat
    fill!(wt, zero(T))
    k0 = argmin(abs.(depth .- depth_shlf))
    wt[k0] = one(T)
    return wt
end

"""
    shelf_weights_interp!(wt, depth, depth_shlf)

Linear interpolation weights between the two layers bracketing `depth_shlf`.
Clamps to the shallowest or deepest layer when out of range.
"""
function shelf_weights_interp!(
    wt         ::AbstractVector{T},
    depth      ::AbstractVector{T},
    depth_shlf ::T,
) where T<:AbstractFloat
    n = length(depth)
    fill!(wt, zero(T))

    if depth_shlf <= depth[1]
        wt[1] = one(T)
        return wt
    elseif depth_shlf >= depth[n]
        wt[n] = one(T)
        return wt
    end

    # Find first index where depth[j] >= depth_shlf
    j = findfirst(>=(depth_shlf), depth)
    if j === nothing
        wt[n] = one(T)
        return wt
    end

    α     = (depth_shlf - depth[j-1]) / (depth[j] - depth[j-1])
    wt[j-1] = one(T) - α
    wt[j]   = α
    return wt
end

# ------------------------------------------------------------
# Shelf-base depth at each (i,j) point
# ------------------------------------------------------------

const ρ_ice_sw_default = 918.0 / 1028.0   # used if no explicit ratio supplied

"""
    shelf_base_depth(H_ice, z_bed, f_grnd, z_sl;
                     depth_mode = :shlf,
                     depth_const = 2000.0,
                     rho_ice_sw  = 918/1028) → depth_2d

Compute the representative depth [m, positive downward] at each grid point.

`depth_mode` options (matching the Fortran `interp_depth` parameter):
- `:shlf`  — floating ice → shelf-base depth; grounded ice → ocean column depth;
             open ocean → `depth_const`.
- `:bed`   — ocean column depth everywhere (`z_sl - z_bed`).
- `:const` — `depth_const` everywhere.
"""
function shelf_base_depth(
    H_ice      ::AbstractMatrix{T},
    z_bed      ::AbstractMatrix{T},
    f_grnd     ::AbstractMatrix{T},
    z_sl       ::AbstractMatrix{T};
    depth_mode ::Symbol = :shlf,
    depth_const::T      = T(2000),
    rho_ice_sw ::T      = T(ρ_ice_sw_default),
) where T<:AbstractFloat
    nx, ny  = size(H_ice)
    out     = similar(H_ice)

    @inbounds for j in 1:ny, i in 1:nx
        if depth_mode === :shlf
            if H_ice[i,j] > zero(T) && f_grnd[i,j] < one(T)
                out[i,j] = H_ice[i,j] * rho_ice_sw
            elseif H_ice[i,j] > zero(T) && f_grnd[i,j] == one(T)
                out[i,j] = z_sl[i,j] - z_bed[i,j]
            else
                out[i,j] = depth_const
            end
        elseif depth_mode === :bed
            out[i,j] = z_sl[i,j] - z_bed[i,j]
        else  # :const
            out[i,j] = depth_const
        end
    end
    return out
end

# ------------------------------------------------------------
# High-level: interpolate a single 3-D field to 2-D
# ------------------------------------------------------------

"""
    interp_to_shelf(var3d, H_ice, z_bed, f_grnd, z_sl, depth;
                    depth_mode    = :shlf,
                    interp_method = :interp,
                    depth_const   = 2000.0,
                    depth_min     = 0.0,
                    depth_max     = 500.0,
                    rho_ice_sw    = 918/1028) → var2d

Vertically interpolate `var3d[nx, ny, nz]` to the shelf-base depth at each
grid point and return the resulting 2-D field.

`interp_method` options (matching the Fortran `interp_method` parameter):
- `:mean`   — equal weighting of all layers in `[depth_min, depth_max]`
- `:layer`  — nearest layer to the target depth
- `:interp` — linear interpolation between the two bracketing layers (default)

`depth_mode` controls how the target depth is determined — see
[`shelf_base_depth`](@ref).
"""
function interp_to_shelf(
    var3d        ::AbstractArray{T,3},
    H_ice        ::AbstractMatrix{T},
    z_bed        ::AbstractMatrix{T},
    f_grnd       ::AbstractMatrix{T},
    z_sl         ::AbstractMatrix{T},
    depth        ::AbstractVector{T};
    depth_mode   ::Symbol = :shlf,
    interp_method::Symbol = :interp,
    depth_const  ::T      = T(2000),
    depth_min    ::T      = T(0),
    depth_max    ::T      = T(500),
    rho_ice_sw   ::T      = T(ρ_ice_sw_default),
) where T<:AbstractFloat

    nx, ny, nz = size(var3d)
    nz == length(depth) ||
        throw(DimensionMismatch("var3d nz=$nz ≠ length(depth)=$(length(depth))"))

    depth_2d = shelf_base_depth(H_ice, z_bed, f_grnd, z_sl;
                                depth_mode, depth_const, rho_ice_sw)
    wt  = zeros(T, nz)
    out = similar(H_ice)

    @inbounds for j in 1:ny, i in 1:nx
        d_shlf = depth_2d[i,j]

        if interp_method === :mean
            shelf_weights_mean!(wt, depth; depth_min, depth_max)
        elseif interp_method === :layer
            shelf_weights_layer!(wt, depth, d_shlf)
        else
            shelf_weights_interp!(wt, depth, d_shlf)
        end

        s = sum(wt)
        if s <= zero(T)
            error("interp_to_shelf: zero weight sum at ($i,$j); check depth array.")
        end
        wt ./= s

        v = zero(T)
        for k in 1:nz
            v += var3d[i,j,k] * wt[k]
        end
        out[i,j] = v
    end

    return out
end

# ------------------------------------------------------------
# Convenience: interpolate T, S, dT together
# ------------------------------------------------------------

"""
    interp_ocean_to_shelf(to_ann, so_ann, dto_ann, H_ice, z_bed, f_grnd, z_sl, depth;
                          kwargs...) → (T_shlf, S_shlf, dT_shlf)

Interpolate ocean temperature `to_ann`, salinity `so_ann`, and temperature
anomaly `dto_ann` from their 3-D representations to the 2-D shelf-base
depth.  All keyword arguments are forwarded to [`interp_to_shelf`](@ref).

Returns a named tuple `(T_shlf, S_shlf, dT_shlf)` ready to be written into
`mshlf.now`.

Pass 2-D arrays for any of the inputs to bypass vertical interpolation for
that variable (a 2-D array is returned unchanged).
"""
function interp_ocean_to_shelf(
    to_ann ::AbstractArray{T},   # 2-D or 3-D
    so_ann ::AbstractArray{T},
    dto_ann::AbstractArray{T},
    H_ice  ::AbstractMatrix{T},
    z_bed  ::AbstractMatrix{T},
    f_grnd ::AbstractMatrix{T},
    z_sl   ::AbstractMatrix{T},
    depth  ::AbstractVector{T};
    kwargs...,
) where T<:AbstractFloat
    _interp(v) = ndims(v) == 3 ?
        interp_to_shelf(v, H_ice, z_bed, f_grnd, z_sl, depth; kwargs...) : v
    return (T_shlf  = _interp(to_ann),
            S_shlf  = _interp(so_ann),
            dT_shlf = _interp(dto_ann))
end
