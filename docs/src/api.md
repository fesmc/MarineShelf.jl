# API Reference

## Types

### Ocean Mask

```@docs
OceanMask
```

### Freezing Point

```@docs
FreezingPointParams
```

### BMB Methods

```@docs
BMBLinear
BMBQuad
BMBQuadNL
BMBAnomaly
BMBPico
```

`BMBMethod{T}` is a convenience type alias for `Union{BMBLinear{T}, BMBQuad{T}, BMBQuadNL{T}, BMBAnomaly{T}, BMBPico{T}}`.

### Thermal-Forcing Methods

```@docs
TFMethod
TFExternal
TFAbsolute
TFAnomaly
```

### Parameter and State Structs

```@docs
MarineShelfParams
MarineShelfState
MShlf
```

---

## Initialisation

```@docs
build_mask_ocn_ref!
```

---

## Update Entry Points

```@docs
update!
update_geometry!
update_tf!
update_bmb!
```

---

## Physics (scalar kernels)

```@docs
freezing_point
bmb_linear
bmb_quad
bmb_quad_nl
bmb_anom
```

---

## Thermal Forcing (in-place)

```@docs
freezing_point!
thermal_forcing_abs!
thermal_forcing_anom!
```

---

## Basin Utilities

```@docs
basin_mean!
basin_mean_icefree!
apply_basin_correction!
```

---

## Limits

```@docs
apply_bmb_limits!
```

---

## Restart I/O

```@docs
write_restart
read_restart!
```

---

## Constants

| Symbol    | Value                        | Units                         | Description                              |
|-----------|------------------------------|-------------------------------|------------------------------------------|
| `ρ_ice`   | 917                          | kg m⁻³                        | Ice density                              |
| `ρ_sw`    | 1028                         | kg m⁻³                        | Sea-water density                        |
| `ρ_w`     | 1000                         | kg m⁻³                        | Fresh-water density                      |
| `cp_o`    | 3974                         | J kg⁻¹ K⁻¹                    | Ocean specific heat capacity             |
| `L_ice`   | 3.34 × 10⁵                  | J kg⁻¹                        | Latent heat of fusion                    |
| `T0`      | 273.15                       | K                             | Celsius to Kelvin offset                 |
| `ρ_ratio` | `ρ_sw / ρ_ice`               | —                             | Density ratio                            |
| `ω`       | `(ρ_sw · cp_o) / (ρ_ice · L_ice)` | K⁻¹                    | Melt efficiency factor                   |
