module MarineShelf

# Julia package for computing the basal mass balance of ice shelves via
# ice–ocean interaction parameterisations.
#
# Design notes:
#   - All 3-D ocean interpolation removed; caller pre-interpolates T/S.
#   - Every numeric struct is parametric in T <: AbstractFloat.
#   - BMB method selection via Julia's type system, not runtime strings.
#   - Thermal-forcing method is explicit: TFAbsolute (default), TFAnomaly
#     (required with BMBAnomaly), or TFExternal (caller pre-sets tf_shlf).
#   - Per-basin corrections pre-built once at construction when basins is
#     supplied, not recomputed each time step.
#
# Quick start:
#   method = BMBQuad{Float32}(γ = Float32(1e4))
#   mshlf  = MShlf(Float32, nx, ny, method;
#                  basins=basins_field, domain="Antarctica",
#                  c_deep=Float32(-50), depth_deep=Float32(2000))
#   build_mask_ocn_ref!(mshlf.now.mask_ocn_ref, regions; domain="Antarctica")
#   # each step:
#   mshlf.now.T_shlf .= T_ocn_at_shelf_base
#   mshlf.now.S_shlf .= S_ocn_at_shelf_base
#   update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx; regions)
#   bmb = mshlf.now.bmb_shlf   # [m/yr], negative = melt

# ------------------------------------------------------------
# Imports
# ------------------------------------------------------------
using NCDatasets

# ------------------------------------------------------------
# Source files (order matters: constants → types → physics → io)
# ------------------------------------------------------------
include("constants.jl")
include("types.jl")
include("geometry.jl")
include("mask.jl")
include("thermal_forcing.jl")
include("bmb.jl")
include("update.jl")
include("io.jl")

# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------

# Types
export OceanMask
export Land, GroundingLine, FloatingLine, Ocean, DeepOcean, Lake
export FreezingPointParams
export BMBLinear, BMBQuad, BMBQuadNL, BMBAnomaly, BMBPico
export TFMethod, TFExternal, TFAbsolute, TFAnomaly
export MarineShelfParams, MarineShelfState, MShlf

# Constants (re-exported for convenience)
export ρ_ice, ρ_sw, ρ_w, cp_o, L_ice, T0, ρ_ratio, ω

# Initialisation helpers
export build_mask_ocn_ref!

# Main update entry points
export update!
export update_geometry!, update_tf!, update_bmb!

# Pure scalar physics (useful for unit tests and offline analysis)
export freezing_point
export bmb_linear, bmb_quad, bmb_quad_nl, bmb_anom

# Thermal forcing helpers (public so callers can set tf_corr externally)
export freezing_point!, thermal_forcing_abs!, thermal_forcing_anom!
export basin_mean!, basin_mean_icefree!
export apply_basin_correction!

# Limits
export apply_bmb_limits!

# I/O
export write_restart, read_restart!

end # module MarineShelf
