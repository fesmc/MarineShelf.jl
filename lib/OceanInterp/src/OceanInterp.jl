"""
    OceanInterp

Pre-processing utilities for preparing ocean fields for use with MarineShelf.jl.

The typical workflow before running the marine-shelf interface model is:

1. **Extrapolate** 3-D ocean fields into ice-shelf cavities (where the ocean
   model has no data):

   ```julia
   ocn_variable_extrapolation!(T_ocn, H_ice, basin_mask, depth, z_bed)
   ocn_variable_extrapolation!(S_ocn, H_ice, basin_mask, depth, z_bed)
   ```

2. **Interpolate** the 3-D fields to the 2-D shelf-base representation that
   MarineShelf expects:

   ```julia
   fields = interp_ocean_to_shelf(T_ocn, S_ocn, dT_ocn,
                                   H_ice, z_bed, f_grnd, z_sl, depth;
                                   depth_mode=:shlf, interp_method=:interp)
   mshlf.now.T_shlf  .= fields.T_shlf
   mshlf.now.S_shlf  .= fields.S_shlf
   mshlf.now.dT_shlf .= fields.dT_shlf
   ```

These functions mirror the public interface of the Fortran `marine_shelf`
module (`ocn_variable_extrapolation`, `marshelf_update_shelf_3D`,
`marshelf_interp_shelf`) but are independent of the time-step loop and
called explicitly by the user.
"""
module OceanInterp

include("extrapolation.jl")
include("vertical_interp.jl")

export ocn_variable_extrapolation!

export shelf_weights_mean!
export shelf_weights_layer!
export shelf_weights_interp!

export shelf_base_depth
export interp_to_shelf
export interp_ocean_to_shelf

end # module OceanInterp
