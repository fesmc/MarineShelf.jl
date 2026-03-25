# ============================================================
#  update.jl
#  Top-level orchestration: update!, update_geometry!,
#  update_tf!, update_bmb!
# ============================================================

# ------------------------------------------------------------
# Full update (geometry + TF + BMB in one call)
# ------------------------------------------------------------

"""
    update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx)

Run a complete marine-shelf update step:

1. Update shelf-base geometry (`z_base`, `slope_base`).
2. Rebuild the ocean mask.
3. Compute the freezing point and thermal forcing.
4. Compute basal mass balance.
5. Apply physical limits.

## Arguments
- `H_ice`   — ice thickness [m]
- `z_bed`   — bed elevation relative to sea level [m]
- `f_grnd`  — grounded fraction [0–1]
- `basins`  — integer basin index per grid point (Float as in Yelmo convention)
- `z_sl`    — sea level [m]
- `dx`      — grid resolution [m]

## Precondition
`mshlf.now.T_shlf`, `mshlf.now.S_shlf`, and (for `BMBAnomaly`)
`mshlf.now.dT_shlf` must be set by the caller before invoking `update!`.
Likewise, `mshlf.now.tf_corr` should be pre-loaded if an external
data correction is desired.
"""
function update!(
    mshlf ::MarineShelf{T},
    H_ice ::AbstractMatrix{T},
    z_bed ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    z_sl  ::AbstractMatrix{T},
    dx    ::T,
) where T<:AbstractFloat
    update_geometry!(mshlf, H_ice, z_bed, f_grnd, z_sl, dx)
    update_tf!(mshlf, H_ice, f_grnd, basins)
    update_bmb!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx)
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
    mshlf ::MarineShelf{T},
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
    update_tf!(mshlf, H_ice, f_grnd, basins)

Compute `now.T_fp_shlf`, `now.tf_corr_basin`, and `now.tf_shlf`.

For `BMBAnomaly` the thermal forcing is assembled from `dT_shlf`;
for all other methods it is assembled from `T_shlf - T_fp`.
"""
function update_tf!(
    mshlf ::MarineShelf{T},
    H_ice ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
) where T<:AbstractFloat
    par = mshlf.par
    now = mshlf.now

    # Apply per-basin TF correction into tf_corr_basin
    fill!(now.tf_corr_basin, zero(T))
    apply_basin_correction!(now.tf_corr_basin, basins, par.basin_tf_corr)

    _update_tf!(now, par.fp, par.bmb_method)
    return mshlf
end

# Dispatch: anomaly method uses dT_shlf directly as TF
function _update_tf!(
    now::MarineShelfState{T},
    fp ::FreezingPointParams{T},
    ::BMBAnomaly,
) where T<:AbstractFloat
    thermal_forcing_anom!(now.tf_shlf, now.dT_shlf,
                          now.tf_corr, now.tf_corr_basin)
end

# Default: all other methods use T_shlf - T_fp
function _update_tf!(
    now::MarineShelfState{T},
    fp ::FreezingPointParams{T},
    ::BMBMethod,
) where T<:AbstractFloat
    freezing_point!(now.T_fp_shlf, now.S_shlf, now.z_base, fp)
    thermal_forcing_abs!(now.tf_shlf, now.T_shlf, now.T_fp_shlf,
                         now.tf_corr, now.tf_corr_basin)
end

# ------------------------------------------------------------
# Step 3: BMB
# ------------------------------------------------------------

"""
    update_bmb!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx)

Compute `now.bmb_shlf` via the configured `BMBMethod`, then apply
per-basin corrections and physical limits.
"""
function update_bmb!(
    mshlf ::MarineShelf{T},
    H_ice ::AbstractMatrix{T},
    z_bed ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    z_sl  ::AbstractMatrix{T},
    dx    ::T,
) where T<:AbstractFloat
    par = mshlf.par
    now = mshlf.now

    # Populate basin-mean fields where the method needs them
    _populate_basin_fields!(now, f_grnd, basins, H_ice, par.bmb_method)

    # Compute raw BMB via dispatched method
    compute_bmb!(now.bmb_shlf, now, par.bmb_method, now.slope_base, dx)

    # Add per-basin BMB correction
    fill!(now.bmb_corr, zero(T))
    apply_basin_correction!(now.bmb_corr, basins, par.basin_bmb_corr)
    @. now.bmb_shlf += now.bmb_corr

    # Physical limits
    apply_bmb_limits!(now.bmb_shlf, H_ice, now.mask_ocn, z_bed, z_sl, par)

    return mshlf
end

# Basin-mean pre-computation — only the methods that need it

function _populate_basin_fields!(
    now   ::MarineShelfState{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    H_ice ::AbstractMatrix{T},
    ::BMBQuadNL{T},
) where T<:AbstractFloat
    # quad-nl needs tf_basin
    basin_mean!(now.tf_basin, now.tf_shlf, f_grnd, basins, H_ice, now.mask_ocn)
end

function _populate_basin_fields!(
    now   ::MarineShelfState{T},
    f_grnd::AbstractMatrix{T},
    basins::AbstractMatrix{T},
    H_ice ::AbstractMatrix{T},
    ::BMBPico{T},
) where T<:AbstractFloat
    # PICO needs ice-free basin means for T and S (far-field forcing)
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
