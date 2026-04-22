# ============================================================
#  types.jl
#  All public and internal types for MarineShelf.
#  Every numeric field is parametric in T <: AbstractFloat so
#  that the working precision is chosen at construction time.
# ============================================================

# ------------------------------------------------------------
# Ocean mask
# ------------------------------------------------------------

"""
Classification of each grid point for marine-shelf calculations.

Values match the legacy integer codes so NetCDF files remain compatible:
- `Land`          = 0
- `GroundingLine` = 1  (grounded side, including partially floating)
- `FloatingLine`  = 2  (first fully-floating ring adjacent to GL)
- `Ocean`         = 3  (open ocean connected to the domain boundary)
- `DeepOcean`     = 4  (open ocean beyond the continental-shelf break)
- `Lake`          = 5  (floating/ocean points isolated from open ocean)
"""
@enum OceanMask::Int8 begin
    Land          = 0
    GroundingLine = 1
    FloatingLine  = 2
    Ocean         = 3
    DeepOcean     = 4
    Lake          = 5
end

# ------------------------------------------------------------
# Freezing-point parameters
# ------------------------------------------------------------

"""
    FreezingPointParams{T}

Linear pressure- and salinity-dependent freezing point (Favier et al. 2019, Eq. 3):

    T_fp = λ1·S + λ2 + λ3·z_base + T0

Default values are standard ISMIP6 coefficients.
"""
Base.@kwdef struct FreezingPointParams{T<:AbstractFloat}
    λ1::T = T(-0.0573)    # °C / PSU
    λ2::T = T( 0.0832)    # °C
    λ3::T = T(-7.61e-4)   # °C / m
end

FreezingPointParams(::Type{T}) where T = FreezingPointParams{T}()

# ------------------------------------------------------------
# BMB method parameter structs  (one per method)
# ------------------------------------------------------------

"""
    BMBLinear{T}

Linear melt parameterisation (Favier et al. 2019, Eq. 2):

    bmb = -γ · ω · TF

If `γ_prime > 0`, γ is scaled by `γ_prime · sin(θ)` (slope-dependent variant).
"""
Base.@kwdef struct BMBLinear{T<:AbstractFloat}
    γ::T
    γ_prime::T = zero(T)
end

"""
    BMBQuad{T}

Local quadratic parameterisation (Favier et al. 2019, Eq. 4):

    bmb = -γ · ω² · TF²

Slope-dependent if `γ_prime > 0`.
"""
Base.@kwdef struct BMBQuad{T<:AbstractFloat}
    γ::T
    γ_prime::T = zero(T)
end

"""
    BMBQuadNL{T}

Non-local (basin-averaged) quadratic parameterisation (Favier et al. 2019, Eq. 5):

    bmb = -γ · ω² · TF · TF_basin

Where `TF_basin` is the basin-mean thermal forcing.
Slope-dependent if `γ_prime > 0`.
"""
Base.@kwdef struct BMBQuadNL{T<:AbstractFloat}
    γ::T
    γ_prime::T = zero(T)
end

"""
    BMBAnomaly{T}

Anomaly-based parameterisation:

    bmb = bmb_ref·c_grz - κ_grz·TF           (grounding zone)
    bmb = above · f_shlf_grz                  (open shelf)

`TF` here is the temperature anomaly `dT_shlf` relative to a reference state,
so this method does NOT require an absolute ocean temperature.  Use with
`tf_method = TFAnomaly()`.
"""
Base.@kwdef struct BMBAnomaly{T<:AbstractFloat}
    κ_grz::T         # m yr⁻¹ K⁻¹ — grounding-zone melt sensitivity
    c_grz::T         # —           — scalar coefficient, grounding zone
    f_grz_shlf::T    # —           — ratio of GL to open-shelf melt rate (≠ 0)
    grz_length::T    # km          — grounding-zone length scale
end

"""
    BMBPico{T}

Placeholder for the PICO box-model method.  PICO state lives in a separate
`Pico` object; the params here are just a tag used for dispatch.
"""
struct BMBPico{T<:AbstractFloat} end

"""Union of all supported BMB method parameter types."""
const BMBMethod{T} = Union{BMBLinear{T}, BMBQuad{T}, BMBQuadNL{T},
                           BMBAnomaly{T}, BMBPico{T}}

# ------------------------------------------------------------
# Thermal-forcing method types  (tf_method)
# ------------------------------------------------------------

"""Abstract supertype for all thermal-forcing computation methods."""
abstract type TFMethod end

"""
    TFExternal

`tf_method = 0`: `tf_shlf` is set by the caller externally before calling
`update!`.  `update_tf!` is a no-op and leaves `tf_shlf` untouched.
"""
struct TFExternal <: TFMethod end

"""
    TFAbsolute

`tf_method = 1`: thermal forcing computed as

    TF = T_shlf - T_fp + tf_corr + tf_corr_basin

Default for all methods except `BMBAnomaly`.
"""
struct TFAbsolute <: TFMethod end

"""
    TFAnomaly

`tf_method = 2`: thermal forcing assembled from the temperature anomaly:

    TF = dT_shlf + tf_corr + tf_corr_basin

Required when `bmb_method = BMBAnomaly`; automatically selected as the
default when that method is chosen.
"""
struct TFAnomaly <: TFMethod end

# Choose sensible defaults: anomaly method requires TFAnomaly
_default_tf_method(::BMBAnomaly) = TFAnomaly()
_default_tf_method(::BMBMethod)  = TFAbsolute()

# ------------------------------------------------------------
# Top-level parameter struct
# ------------------------------------------------------------

"""
    MarineShelfParams{T, M<:BMBMethod{T}, F<:TFMethod}

All parameters needed to run `update!`.

# Fields
- `bmb_method`      — one of `BMBLinear`, `BMBQuad`, `BMBQuadNL`,
                       `BMBAnomaly`, `BMBPico`
- `tf_method`       — one of `TFExternal`, `TFAbsolute`, `TFAnomaly`
- `fp`              — freezing-point coefficients
- `c_deep`          — melt rate applied to deep-ocean points [m/yr]
- `depth_deep`      — ocean depth threshold for `c_deep` [m]
- `bmb_max`         — maximum allowed refreezing rate [m/yr]
- `find_ocean`      — flood-fill to distinguish open ocean from lakes
- `domain`          — model domain string (e.g. `"Antarctica"`, `"Greenland"`)
- `basin_tf_corr`   — per-basin additive correction to thermal forcing
- `basin_bmb_corr`  — per-basin additive correction to bmb
"""
struct MarineShelfParams{T<:AbstractFloat, M<:BMBMethod, F<:TFMethod}
    bmb_method::M
    tf_method::F
    fp::FreezingPointParams{T}
    c_deep::T
    depth_deep::T
    bmb_max::T
    find_ocean::Bool
    domain::String
    basin_tf_corr::Dict{Int,T}
    basin_bmb_corr::Dict{Int,T}
end

"""
    MarineShelfParams(T, bmb_method; kwargs...)

Convenience constructor.  `T` sets the working precision.
"""
function MarineShelfParams(
    ::Type{T},
    bmb_method::M;
    tf_method      = _default_tf_method(bmb_method),
    fp             = FreezingPointParams(T),
    c_deep         = T(-50),
    depth_deep     = T(2000),
    bmb_max        = T(0),
    find_ocean     = true,
    domain         = "",
    basin_tf_corr  = Dict{Int,T}(),
    basin_bmb_corr = Dict{Int,T}(),
) where {T<:AbstractFloat, M<:BMBMethod}
    if tf_method isa TFAbsolute && bmb_method isa BMBAnomaly
        error("tf_method=TFAbsolute() is incompatible with bmb_method=BMBAnomaly. " *
              "Use TFAnomaly() or TFExternal().")
    end
    F = typeof(tf_method)
    return MarineShelfParams{T,M,F}(
        bmb_method, tf_method, fp,
        T(c_deep), T(depth_deep), T(bmb_max),
        find_ocean, domain,
        basin_tf_corr, basin_bmb_corr,
    )
end

# ------------------------------------------------------------
# State struct
# ------------------------------------------------------------

"""
    MarineShelfState{T, A<:AbstractMatrix{T}}

All 2-D fields maintained across time steps.

## Fields set externally before calling `update!`
- `T_shlf`   — ocean temperature at the shelf base [K]
- `S_shlf`   — ocean salinity at the shelf base [PSU]
- `dT_shlf`  — temperature anomaly relative to a reference state [K]
- `dS_shlf`  — salinity anomaly relative to a reference state [PSU]

## Fields updated internally by `update!`
- `z_base`, `slope_base` — shelf-base geometry
- `T_fp_shlf`            — pressure-dependent freezing point [K]
- `tf_shlf`              — thermal forcing [K]
- `tf_basin`             — basin-mean thermal forcing [K]  (quad-nl / pico)
- `tf_corr`              — external data correction to TF [K]  (e.g. ISMIP6)
- `tf_corr_basin`        — per-basin index correction to TF [K]  (set once at init)
- `T_shlf_basin`, `S_shlf_basin` — basin-mean T/S (pico)
- `bmb_shlf`             — basal mass balance [m/yr], negative = melt
- `bmb_ref`              — reference BMB for the anomaly method [m/yr]
- `bmb_corr`             — per-basin additive correction to bmb [m/yr]  (set once at init)
- `mask_ocn_ref`         — static reference ocean mask (set at init)
- `mask_ocn`             — dynamic ocean mask (updated each step)
"""
mutable struct MarineShelfState{T<:AbstractFloat, A<:AbstractMatrix{T}}
    # --- Inputs (caller-owned, written before update!) -----------------
    T_shlf::A
    S_shlf::A
    dT_shlf::A
    dS_shlf::A

    # --- Geometry (computed by update_geometry!) -----------------------
    z_base::A
    slope_base::A

    # --- Thermal forcing (computed by update_tf!) ----------------------
    T_fp_shlf::A
    tf_shlf::A
    tf_basin::A
    tf_corr::A
    tf_corr_basin::A

    # --- Basin-mean fields for PICO / quad-nl --------------------------
    T_shlf_basin::A
    S_shlf_basin::A

    # --- BMB output (computed by update_bmb!) --------------------------
    bmb_shlf::A
    bmb_ref::A
    bmb_corr::A

    # --- Masks ---------------------------------------------------------
    mask_ocn_ref::Matrix{OceanMask}
    mask_ocn::Matrix{OceanMask}
end

"""
    MarineShelfState(T, nx, ny; array_constructor=zeros)

Allocate a zero-filled state on an `nx × ny` grid.
`array_constructor` can be replaced with e.g. `CUDA.zeros` for GPU arrays.
"""
function MarineShelfState(
    ::Type{T}, nx::Int, ny::Int;
    array_constructor = (T, nx, ny) -> zeros(T, nx, ny),
) where T<:AbstractFloat
    A = array_constructor(T, nx, ny)
    newfield() = array_constructor(T, nx, ny)
    return MarineShelfState{T, typeof(A)}(
        newfield(), newfield(), newfield(), newfield(),   # T_shlf, S_shlf, dT_shlf, dS_shlf
        newfield(), newfield(),                           # z_base, slope_base
        newfield(), newfield(), newfield(),               # T_fp_shlf, tf_shlf, tf_basin
        newfield(), newfield(),                           # tf_corr, tf_corr_basin
        newfield(), newfield(),                           # T_shlf_basin, S_shlf_basin
        newfield(), newfield(), newfield(),               # bmb_shlf, bmb_ref, bmb_corr
        fill(Ocean, nx, ny),                             # mask_ocn_ref
        fill(Ocean, nx, ny),                             # mask_ocn
    )
end

# ------------------------------------------------------------
# Top-level object
# ------------------------------------------------------------

"""
    MShlf{T, M, A}

Top-level object combining parameters (`par`) and time-evolving state (`now`).
`T` is the working precision, `M` the BMB method type, `A` the array type.

Construct via `MShlf(T, nx, ny, bmb_method; kwargs...)`.
"""
struct MShlf{T<:AbstractFloat, M<:BMBMethod, A<:AbstractMatrix{T}}
    par::MarineShelfParams{T,M}
    now::MarineShelfState{T,A}
end

"""
    MShlf(T, nx, ny, bmb_method; basins=nothing, kwargs...)

Allocate a marine-shelf model on an `nx × ny` grid at precision `T`.

If per-basin corrections (`basin_tf_corr` or `basin_bmb_corr`) are
specified, `basins` must be provided so the correction fields can be
pre-built once at construction (matching the Fortran `marshelf_init`
behaviour).  The `tf_corr_basin` and `bmb_corr` fields in state are
populated once here and are not modified during the time loop.

`kwargs` are forwarded to `MarineShelfParams`.
"""
function MShlf(
    ::Type{T}, nx::Int, ny::Int, bmb_method::M;
    basins            = nothing,
    array_constructor = (T, nx, ny) -> zeros(T, nx, ny),
    kwargs...
) where {T<:AbstractFloat, M<:BMBMethod}
    par = MarineShelfParams(T, bmb_method; kwargs...)
    now = MarineShelfState(T, nx, ny; array_constructor)

    has_corr = !isempty(par.basin_tf_corr) || !isempty(par.basin_bmb_corr)
    if has_corr && isnothing(basins)
        error("MShlf: `basins` must be provided when " *
              "`basin_tf_corr` or `basin_bmb_corr` are non-empty.")
    end
    if !isnothing(basins)
        apply_basin_correction!(now.tf_corr_basin, basins, par.basin_tf_corr)
        apply_basin_correction!(now.bmb_corr,      basins, par.basin_bmb_corr)
    end

    return MShlf{T, M, typeof(now.bmb_shlf)}(par, now)
end
