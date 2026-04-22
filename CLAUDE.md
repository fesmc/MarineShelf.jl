# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Instantiate dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run tests directly
julia --project=. test/runtests.jl

# Load the package in a REPL
julia --project=. -e 'using MarineShelf'
```

## Architecture

MarineShelf.jl computes the basal mass balance (BMB) of ice shelves via ice–ocean interaction parameterisations. It is designed as a component driven by an ice-sheet model; the caller pre-interpolates ocean T/S to the shelf-base depth before passing them in.

**Source file load order** (defined in `src/MarineShelf.jl`) must be respected:
`constants.jl` → `types.jl` → `geometry.jl` → `mask.jl` → `thermal_forcing.jl` → `bmb.jl` → `update.jl` → `io.jl`

### Type hierarchy

- `MShlf{T,M,A}` — top-level object: holds `par::MarineShelfParams` and `now::MarineShelfState`. Named `MShlf` (not `MarineShelf`) to avoid a name clash with the enclosing module when using `using MarineShelf`.
- `MarineShelfParams{T,M,F}` — static config: `bmb_method`, `tf_method`, freezing-point coefficients, depth limits, `domain`, per-basin correction dicts
- `MarineShelfState{T,A}` — all mutable 2-D fields (inputs, geometry, TF, BMB, masks)
- `BMBMethod{T}` — union of `BMBLinear`, `BMBQuad`, `BMBQuadNL`, `BMBAnomaly`, `BMBPico`; method selection entirely via Julia dispatch
- `TFMethod` — abstract type with three concrete singletons:
  - `TFExternal()` — caller pre-sets `tf_shlf`; `update_tf!` is a no-op
  - `TFAbsolute()` — `TF = T_shlf - T_fp + corrections` (default for all BMB methods except `BMBAnomaly`)
  - `TFAnomaly()` — `TF = dT_shlf + corrections` (default for `BMBAnomaly`, required)
- `OceanMask` — `Int8` enum: `Land`(0), `GroundingLine`(1), `FloatingLine`(2), `Ocean`(3), `DeepOcean`(4), `Lake`(5)

Every numeric struct is parametric in `T <: AbstractFloat`.

### Update pipeline (`update.jl`)

`update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx; regions=nothing)` runs three sequential steps:

1. **`update_geometry!`** — `z_base`, `slope_base`, and `mask_ocn` via flood-fill from `mask_ocn_ref`
2. **`update_tf!`** — dispatches on `par.tf_method`; `TFExternal` is a no-op, `TFAbsolute` computes freezing point then TF, `TFAnomaly` uses `dT_shlf`
3. **`update_bmb!`** — pre-populates basin-mean fields if needed, calls dispatched `compute_bmb!`, optionally applies PICO/Antarctica blend (when `par.domain == "Antarctica"` and `regions` is supplied), adds pre-built `bmb_corr`, applies physical limits

### Basin corrections (static, set at construction)

`bmb_corr` and `tf_corr_basin` in state are **pre-built once** when `MShlf` is constructed (matching the Fortran `marshelf_init` behaviour). Pass `basins` to the constructor along with `basin_bmb_corr` / `basin_tf_corr` dicts. These fields are not touched during the time loop. Providing non-empty correction dicts without `basins` raises an error.

### Caller responsibilities (preconditions before `update!`)

- `mshlf.now.T_shlf`, `mshlf.now.S_shlf` — ocean T [K] and S [PSU] pre-interpolated to shelf-base depth
- `mshlf.now.dT_shlf` — temperature anomaly [K] (only for `TFAnomaly`)
- `mshlf.now.dS_shlf` — salinity anomaly [PSU] (diagnostic, not used in physics)
- `mshlf.now.tf_corr` — external TF correction (e.g. ISMIP6), if desired
- `mshlf.now.tf_shlf` — when `tf_method = TFExternal()`, must be set by caller
- `mshlf.now.mask_ocn_ref` — set once at init via `build_mask_ocn_ref!(mshlf.now.mask_ocn_ref, regions; domain="Antarctica")`

### BMB parameterisations (`bmb.jl`)

| Type | Formula |
|------|---------|
| `BMBLinear` | `−γ·ω·TF` |
| `BMBQuad` | `−γ·ω²·TF²` |
| `BMBQuadNL` | `−γ·ω²·TF·TF_basin` (non-local; sign flip when both TF negative) |
| `BMBAnomaly` | grounding-zone/open-shelf blend using `dT_shlf` as TF |
| `BMBPico` | no-op; BMB populated externally by a PICO sub-model |

All methods support an optional slope-dependent exchange velocity when `γ_prime > 0`.

### I/O (`io.jl`)

`write_restart` / `read_restart!` use NCDatasets.jl. Parameters are not stored — only state fields including `dS_shlf`. Coordinates stored in km; data fields stored as `Float32`.

## OceanInterp subpackage (`lib/OceanInterp/`)

Standalone utilities for preparing 3-D ocean data before driving MarineShelf. Typical workflow:

1. **Extrapolate** into ice-shelf cavities: `ocn_variable_extrapolation!(var3d, H_ice, basin_mask, depth, z_bed)`
2. **Interpolate** to 2-D shelf representation: `interp_ocean_to_shelf(T, S, dT, H_ice, z_bed, f_grnd, z_sl, depth; depth_mode, interp_method, ...)`

Three vertical interpolation methods (`interp_method`): `:mean` (depth-range average), `:layer` (nearest layer), `:interp` (linear; default). Three depth selection modes (`depth_mode`): `:shlf` (shelf-base depth; default), `:bed` (ocean column depth), `:const` (fixed depth).

These functions are called explicitly by the caller; they are not invoked by `update!`.
