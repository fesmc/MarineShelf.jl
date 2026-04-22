# OceanInterp

The `OceanInterp` subpackage provides pre-processing utilities for preparing ocean fields before passing them to MarineShelf.jl. It mirrors the public interface of the Fortran `marine_shelf` module's interpolation routines (`ocn_variable_extrapolation`, `marshelf_interp_shelf`) but is called explicitly by the user, separate from the time-step loop.

## Installation

`OceanInterp` is a standalone package in the `lib/` subdirectory:

```julia
using Pkg
Pkg.develop(path = "lib/OceanInterp")
using OceanInterp
```

## Typical Workflow

```julia
using OceanInterp, MarineShelf

# Step 1 — extrapolate 3-D fields into ice-shelf cavities
ocn_variable_extrapolation!(T_ocn, H_ice, basin_mask, depth, z_bed)
ocn_variable_extrapolation!(S_ocn, H_ice, basin_mask, depth, z_bed)

# Step 2 — interpolate to the 2-D shelf-base representation
fields = interp_ocean_to_shelf(T_ocn, S_ocn, dT_ocn,
                                H_ice, z_bed, f_grnd, z_sl, depth;
                                depth_mode    = :shlf,
                                interp_method = :interp)

mshlf.now.T_shlf  .= fields.T_shlf
mshlf.now.S_shlf  .= fields.S_shlf
mshlf.now.dT_shlf .= fields.dT_shlf
```

---

## API

### Cavity Extrapolation

```@docs
OceanInterp.ocn_variable_extrapolation!
```

### Depth Selection

```@docs
OceanInterp.shelf_base_depth
```

`depth_mode` controls how the representative depth is determined at each grid point:

| Mode     | Description                                                                     |
|----------|---------------------------------------------------------------------------------|
| `:shlf`  | Floating ice → shelf-base depth; grounded ice → ocean column; open ocean → `depth_const` |
| `:bed`   | Ocean column depth (`z_sl - z_bed`) everywhere                                  |
| `:const` | `depth_const` everywhere                                                        |

### Vertical Weight Functions

```@docs
OceanInterp.shelf_weights_mean!
OceanInterp.shelf_weights_layer!
OceanInterp.shelf_weights_interp!
```

`interp_method` options passed to `interp_to_shelf` / `interp_ocean_to_shelf`:

| Method    | Description                                                   |
|-----------|---------------------------------------------------------------|
| `:mean`   | Equal weighting of all layers in `[depth_min, depth_max]`     |
| `:layer`  | All weight on the single nearest layer                        |
| `:interp` | Linear interpolation between the two bracketing layers (default) |

### Field Interpolation

```@docs
OceanInterp.interp_to_shelf
OceanInterp.interp_ocean_to_shelf
```
