"""
    MarineShelf

Julia package for computing the basal mass balance of ice shelves via
ice–ocean interaction parameterisations.

## Design notes

- All 3-D ocean interpolation has been removed.  The caller is responsible
  for providing `T_shlf` and `S_shlf` pre-interpolated to the relevant
  depth of the shelf base.
- Every numeric struct is parametric in `T <: AbstractFloat`, so the
  working precision is chosen at construction time (e.g. `Float32` for
  performance, `Float64` for reference runs).
- BMB method selection is done via Julia's type system; there is no
  runtime string dispatch.

## Quick start

```julia
using MarineShelf

# Choose method and precision
method = BMBQuad{Float32}(γ = 1e4f0)
mshlf  = MarineShelf(Float32, nx, ny, method;
                     c_deep=Float32(-50), depth_deep=Float32(2000),
                     bmb_max=Float32(0))

# Initialise the reference ocean mask (once, before the time loop)
build_mask_ocn_ref!(mshlf.now.mask_ocn_ref, regions; domain="Antarctica")

# Each time step: set inputs, then update
mshlf.now.T_shlf  .= T_ocn_at_shelf_base   # pre-interpolated [K]
mshlf.now.S_shlf  .= S_ocn_at_shelf_base   # pre-interpolated [PSU]
update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx)

# Read result
bmb = mshlf.now.bmb_shlf   # [m/yr], negative = melt
```
"""
module MarineShelf

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
export MarineShelfParams, MarineShelfState, MarineShelf

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
