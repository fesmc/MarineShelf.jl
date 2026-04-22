# User Guide

## Workflow

A typical integration looks like this:

```
Allocate MShlf (once)
    └── MarineShelfParams  ← BMB method, TF method, basin corrections, limits
    └── MarineShelfState   ← 2-D fields (T_shlf, S_shlf, z_base, bmb_shlf, …)

Build mask_ocn_ref (once, after allocating)

Per time step:
    Set now.T_shlf, now.S_shlf  (from your ocean driver)
    update!(mshlf, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    read mshlf.now.bmb_shlf
```

`update!` performs three sub-steps in order:

1. **Geometry** — shelf-base elevation `z_base`, local slope `slope_base`, and the dynamic ocean mask `mask_ocn`.
2. **Thermal forcing** — freezing point `T_fp_shlf`, then `tf_shlf = T_shlf - T_fp + corrections`.
3. **BMB** — raw parameterisation, optional PICO/Antarctica blend, per-basin `bmb_corr`, physical limits.

You can call each sub-step independently via `update_geometry!`, `update_tf!`, `update_bmb!`.

---

## BMB Methods

All methods require a parameter struct parametric in your working type `T`:

### `BMBLinear{T}` — linear parameterisation

```
bmb = -γ · ω · TF
```

```julia
method = BMBLinear{Float64}(γ = 1e4)
```

### `BMBQuad{T}` — local quadratic (default choice)

```
bmb = -γ · ω² · TF²
```

```julia
method = BMBQuad{Float64}(γ = 1e4)
```

### `BMBQuadNL{T}` — non-local (basin-averaged) quadratic

```
bmb = -γ · ω² · TF · TF_basin
```
where `TF_basin` is the mean thermal forcing over all floating-ice points in each basin. When both `TF` and `TF_basin` are negative the sign is flipped to allow accretion.

```julia
method = BMBQuadNL{Float64}(γ = 1e4)
```

### `BMBAnomaly{T}` — anomaly-based

```
bmb_gl    = c_grz · bmb_ref - κ_grz · dT_shlf
bmb_shelf = bmb_gl / f_grz_shlf
```
At grounding-line / floating-line points a resolution-weighted blend between the two is returned. Requires `tf_method = TFAnomaly()` (selected automatically); `now.bmb_ref` must be loaded by the caller.

```julia
method = BMBAnomaly{Float64}(
    κ_grz      = 8.0,     # m yr⁻¹ K⁻¹
    c_grz      = 1.0,
    f_grz_shlf = 4.0,
    grz_length = 20.0,    # km
)
```

### `BMBPico{T}` — PICO placeholder

The full PICO box model lives in a separate object; this struct is a dispatch tag only. BMB must be computed externally and written into `now.bmb_shlf` before `update_bmb!` is called (or `update!` is called). `update_bmb!` then applies limits and per-basin corrections.

```julia
method = BMBPico{Float64}()
```

### Slope-dependent exchange velocity

`BMBLinear`, `BMBQuad`, and `BMBQuadNL` all accept an optional `γ_prime` parameter. When `γ_prime > 0` the exchange velocity is scaled by `γ_prime · sin(θ)` at each grid point.

---

## Thermal-Forcing Methods

### `TFAbsolute()` (default)

```
TF = T_shlf - T_fp(S_shlf, z_base) + tf_corr + tf_corr_basin
```
Use with `BMBLinear`, `BMBQuad`, `BMBQuadNL`, `BMBPico`.

### `TFAnomaly()` (required for `BMBAnomaly`)

```
TF = dT_shlf + tf_corr + tf_corr_basin
```
Set `now.dT_shlf` from your ocean anomaly field.

### `TFExternal()`

`update_tf!` is a no-op. Set `now.tf_shlf` directly before calling `update!`.

---

## Per-Basin Corrections

Pass correction dictionaries at construction; the 2-D correction fields are built once and never recomputed:

```julia
mshlf = MShlf(Float64, nx, ny, method;
              basins         = basins_field,
              basin_tf_corr  = Dict(1 => 0.5, 2 => -0.3),   # TF correction per basin
              basin_bmb_corr = Dict(3 => -2.0))               # BMB correction per basin
```

The `basins` field must be provided whenever either dict is non-empty.

An external 2-D correction field (e.g. loaded from ISMIP6 data) can be added to `now.tf_corr` directly at any time; it is added on top of the basin correction.

---

## Ocean Mask

### Static reference mask

`mask_ocn_ref` is set once after allocation and never touched by `update!`:

```julia
build_mask_ocn_ref!(mshlf.now.mask_ocn_ref, regions; domain = "Antarctica")
# or domain = "Greenland", or omit for generic
```

You can modify `mask_ocn_ref` directly afterwards for custom adjustments.

### Dynamic mask

`mask_ocn` is rebuilt each call to `update_geometry!`. It classifies:

| Value          | Meaning                                             |
|----------------|-----------------------------------------------------|
| `Land`         | fully grounded, no floating neighbours              |
| `GroundingLine`| grounded cell adjacent to floating ice             |
| `FloatingLine` | first fully-floating ring adjacent to grounded ice  |
| `Ocean`        | open ocean connected to domain boundary             |
| `DeepOcean`    | open ocean beyond the continental-shelf break       |
| `Lake`         | ocean/floating points isolated from open ocean      |

The flood-fill that distinguishes `Ocean` from `Lake` can be disabled with `find_ocean=false`.

---

## Physical Limits

Applied at the end of each `update_bmb!`:

1. **Refreezing cap** — `bmb[i,j] = min(bmb[i,j], par.bmb_max)` (default `bmb_max = 0`)
2. **No accretion without ice** — zero out positive BMB where `H_ice == 0`
3. **Deep-ocean forcing** — apply `c_deep` at `DeepOcean` points deeper than `depth_deep`

---

## Restart I/O

```julia
# Write a restart snapshot
write_restart("restart.nc", mshlf, xc, yc, time)

# Append a second snapshot
write_restart("restart.nc", mshlf, xc, yc, time + dt)

# Read back (last slice)
read_restart!(mshlf, "restart.nc"; time_index=nothing)
```

The file stores all state fields needed to warm-restart a run. Parameters are not stored — reload them from your namelist or configuration.
