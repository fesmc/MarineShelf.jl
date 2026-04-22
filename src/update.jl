# ============================================================
#  update.jl
#  Top-level orchestration: update!, update_geometry!,
#  update_tf!, update_bmb!
# ============================================================

# ------------------------------------------------------------
# Full update (geometry + TF + BMB in one call)
# ------------------------------------------------------------

"""
    update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx; regions=nothing)

Run a complete marine-shelf update step:

1. Update shelf-base geometry (`z_base`, `slope_base`).
2. Rebuild the ocean mask.
3. Compute the freezing point and thermal forcing (unless `tf_method = TFExternal()`).
4. Compute basal mass balance.
5. Apply physical limits.

## Arguments
- `H_ice`   — ice thickness [m]
- `z_bed`   — bed elevation relative to sea level [m]
- `f_grnd`  — grounded fraction [0–1]
- `basins`  — integer basin index per grid point (Float as in Yelmo convention)
- `z_sl`    — sea level [m]
- `dx`      — grid resolution [m]

## Keyword arguments
- `regions` — (optional) region mask used for the PICO/Antarctica special case

## Preconditions
`mshlf.now.T_shlf`, `mshlf.now.S_shlf`, and (for `BMBAnomaly`)
`mshlf.now.dT_shlf` must be set by the caller.
`mshlf.now.tf_corr` should be pre-loaded if an external data correction is desired.
When `tf_method = TFExternal()`, `mshlf.now.tf_shlf` must be set by the caller.
"""
function update!(
    mshlf ::MShlf{T},
    H_ice ::AbstractMatrix{T},
    z_bed ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    z_sl  ::AbstractMatrix{T},
    dx    ::T;
    regions::Union{AbstractMatrix,Nothing} = nothing,
) where T<:AbstractFloat
    update_geometry!(mshlf, H_ice, z_bed, f_grnd, z_sl, dx)
    update_tf!(mshlf)
    update_bmb!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx; regions)
    return mshlf
end

# ------------------------------------------------------------
# Step 1: geometry
# ------------------------------------------------------------

"""
    update_geometry!(mshlf, H_ice, z_bed, f_grnd, z_sl, dx)

Update `now.z_base`, `now.slope_base`, and `now.mask_ocn`.
"""
function update_geometry!(
    mshlf ::MShlf{T},
    H_ice ::AbstractMatrix{T},
    z_bed ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    z_sl  ::AbstractMatrix{T},
    dx    ::T,
) where T<:AbstractFloat
    now = mshlf.now
    shelf_base_elevation!(now.z_base, H_ice, f_grnd, z_bed, z_sl)
    shelf_base_slope!(now.slope_base, now.z_base, dx)
    ocean_mask!(now.mask_ocn, f_grnd, now.mask_ocn_ref;
                find_ocean=mshlf.par.find_ocean)
    return mshlf
end

# ------------------------------------------------------------
# Step 2: thermal forcing
# ------------------------------------------------------------

"""
    update_tf!(mshlf)

Compute thermal forcing according to `mshlf.par.tf_method`:

- `TFExternal()` — no-op; caller is responsible for setting `now.tf_shlf`.
- `TFAbsolute()` — computes `T_fp_shlf` then `tf_shlf = T_shlf - T_fp + corrections`.
- `TFAnomaly()`  — computes `tf_shlf = dT_shlf + corrections`.

Per-basin corrections (`tf_corr_basin`) are pre-built at construction and
added here without recomputation.
"""
function update_tf!(mshlf::MShlf)
    _update_tf!(mshlf.now, mshlf.par.fp, mshlf.par.tf_method)
    return mshlf
end

# tf_method = 0: externally supplied — do nothing
function _update_tf!(::MarineShelfState, ::FreezingPointParams, ::TFExternal) end

# tf_method = 1: absolute temperature method
function _update_tf!(
    now::MarineShelfState{T},
    fp ::FreezingPointParams{T},
    ::TFAbsolute,
) where T<:AbstractFloat
    freezing_point!(now.T_fp_shlf, now.S_shlf, now.z_base, fp)
    thermal_forcing_abs!(now.tf_shlf, now.T_shlf, now.T_fp_shlf,
                         now.tf_corr, now.tf_corr_basin)
end

# tf_method = 2: anomaly method
function _update_tf!(
    now::MarineShelfState{T},
    ::FreezingPointParams{T},
    ::TFAnomaly,
) where T<:AbstractFloat
    thermal_forcing_anom!(now.tf_shlf, now.dT_shlf,
                          now.tf_corr, now.tf_corr_basin)
end

# ------------------------------------------------------------
# Step 3: BMB
# ------------------------------------------------------------

"""
    update_bmb!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx; regions=nothing)

Compute `now.bmb_shlf` via the configured `BMBMethod`, then add the
pre-built per-basin correction and apply physical limits.

For `bmb_method = BMBPico()` on the `"Antarctica"` domain, the optional
`regions` matrix is used to blend shelf-margin melt with `c_deep` at
region 2.11 boundaries (matching the Fortran behaviour).
"""
function update_bmb!(
    mshlf ::MShlf{T},
    H_ice ::AbstractMatrix{T},
    z_bed ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    z_sl  ::AbstractMatrix{T},
    dx    ::T;
    regions::Union{AbstractMatrix,Nothing} = nothing,
) where T<:AbstractFloat
    par = mshlf.par
    now = mshlf.now

    # Populate basin-mean fields where the method needs them
    _populate_basin_fields!(now, f_grnd, basins, H_ice, par.bmb_method)

    # Compute raw BMB via dispatched method
    compute_bmb!(now.bmb_shlf, now, par.bmb_method, now.slope_base, dx)

    # PICO / Antarctica special case (Fortran: blend with c_deep at region 2.1)
    if par.bmb_method isa BMBPico && par.domain == "Antarctica" && !isnothing(regions)
        _apply_pico_antarctica!(now.bmb_shlf, regions, par.c_deep)
    end

    # Add pre-built per-basin BMB correction (set once at construction)
    @. now.bmb_shlf += now.bmb_corr

    # Physical limits
    apply_bmb_limits!(now.bmb_shlf, H_ice, now.mask_ocn, z_bed, z_sl, par)

    return mshlf
end

@inline function _apply_pico_antarctica!(
    bmb    ::AbstractMatrix{T},
    regions::AbstractMatrix,
    c_deep ::T,
) where T<:AbstractFloat
    half = T(0.5)
    @inbounds for j in axes(bmb, 2), i in axes(bmb, 1)
        if regions[i,j] == 2.1 || regions[i,j] == T(2.1)
            bmb[i,j] = half * bmb[i,j] + half * c_deep
        end
    end
end

# Basin-mean pre-computation — only the methods that need it

function _populate_basin_fields!(
    now   ::MarineShelfState{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    H_ice ::AbstractMatrix{T},
    ::BMBQuadNL{T},
) where T<:AbstractFloat
    basin_mean!(now.tf_basin, now.tf_shlf, f_grnd, basins, H_ice, now.mask_ocn)
end

function _populate_basin_fields!(
    now   ::MarineShelfState{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    H_ice ::AbstractMatrix{T},
    ::BMBPico{T},
) where T<:AbstractFloat
    basin_mean_icefree!(now.T_shlf_basin, now.T_shlf,
                        f_grnd, basins, H_ice, now.mask_ocn)
    basin_mean_icefree!(now.S_shlf_basin, now.S_shlf,
                        f_grnd, basins, H_ice, now.mask_ocn)
end

# All other methods: nothing to pre-compute
function _populate_basin_fields!(
    ::MarineShelfState{T},
    ::AbstractMatrix{T},
    ::AbstractMatrix{T},
    ::AbstractMatrix{T},
    ::BMBMethod,
) where T<:AbstractFloat
end
