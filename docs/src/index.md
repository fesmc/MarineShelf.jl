# MarineShelf.jl

**MarineShelf.jl** is a Julia package for computing the basal mass balance (BMB) of marine ice shelves via ice–ocean interaction parameterisations. It is a clean Julia port of the `marine_shelf` Fortran module from the [Yelmo](https://github.com/palma-ice/yelmox) ice-sheet model, preserving the same physics and public interface while taking advantage of Julia's type system for method dispatch.

## Features

- Five BMB parameterisations: linear, local quadratic, non-local quadratic, anomaly-based, and PICO placeholder
- Three thermal-forcing methods: absolute temperature, anomaly, or externally supplied
- Pressure-dependent freezing point (ISMIP6 coefficients)
- Per-basin thermal-forcing and BMB corrections, pre-built once at construction
- Ocean mask with flood-fill to classify lakes, grounding-line zones, and deep ocean
- NetCDF restart I/O via [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl)
- Companion [OceanInterp](oceaninterp.md) subpackage for 3-D cavity extrapolation and shelf interpolation

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/fesmc/MarineShelf.jl")
```

## Quick Start

```julia
using MarineShelf

# 1. Choose a BMB method
method = BMBQuad{Float32}(γ = Float32(1e4))

# 2. Allocate the model
mshlf = MShlf(Float32, nx, ny, method;
              domain     = "Antarctica",
              c_deep     = Float32(-50),
              depth_deep = Float32(2000))

# 3. Build the static reference ocean mask (once)
build_mask_ocn_ref!(mshlf.now.mask_ocn_ref, regions; domain = "Antarctica")

# 4. Time loop
for t in timesteps
    # Set ocean fields from your ocean model / interpolation
    mshlf.now.T_shlf .= T_ocn_at_shelf_base
    mshlf.now.S_shlf .= S_ocn_at_shelf_base

    # Run one update step
    update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx)

    # Basal mass balance [m yr⁻¹], negative = melt
    bmb = mshlf.now.bmb_shlf
end
```

See the [User Guide](guide.md) for a complete description of all methods and options.
