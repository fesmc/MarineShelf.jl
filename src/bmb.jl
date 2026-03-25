# ============================================================
#  bmb.jl
#  Basal mass balance parameterisations.
#
#  Each method is dispatched via Julia's type system on the
#  parameter struct rather than a runtime string switch.
#  Pure scalar functions are provided alongside in-place 2-D
#  wrappers so that individual methods are easily unit-tested.
# ============================================================

# ============================================================
# Scalar kernel functions (pure, no allocation)
# ============================================================

"""
    bmb_linear(tf, γ, ω) → T

Linear parameterisation (Favier et al. 2019, Eq. 2):

    bmb = -γ · ω · TF

- `tf` — thermal forcing [K]
- `γ`  — heat exchange velocity [m yr⁻¹]
- `ω`  — `(ρ_sw · cp_o) / (ρ_ice · L_ice)` [K⁻¹]
"""
@inline bmb_linear(tf::T, γ::T, ω_val::T) where T<:AbstractFloat =
    -γ * ω_val * tf

"""
    bmb_quad(tf, γ, ω) → T

Local quadratic parameterisation (Favier et al. 2019, Eq. 4):

    bmb = -γ · ω² · TF²
"""
@inline bmb_quad(tf::T, γ::T, ω_val::T) where T<:AbstractFloat =
    -γ * ω_val^2 * tf^2

"""
    bmb_quad_nl(tf, tf_basin, γ, ω) → T

Non-local quadratic parameterisation (Favier et al. 2019, Eq. 5):

    bmb = -γ · ω² · TF · TF_basin
"""
@inline bmb_quad_nl(tf::T, tf_basin::T, γ::T, ω_val::T) where T<:AbstractFloat =
    -γ * ω_val^2 * tf * tf_basin

"""
    bmb_anom(tf, bmb_ref, mask, p::BMBAnomaly, dx) → T

Anomaly-based parameterisation.

Grounding-zone value:
    bmb_gl = c_grz · bmb_ref - κ_grz · TF

Open-shelf value:
    bmb_shlf = bmb_gl / f_grz_shlf

At grounding-line / floating-line points a resolution-weighted blend
between the two is returned.
"""
@inline function bmb_anom(
    tf      ::T,
    bmb_ref ::T,
    mask    ::OceanMask,
    p       ::BMBAnomaly{T},
    dx      ::T,
) where T<:AbstractFloat
    f_shlf = one(T) / p.f_grz_shlf
    grz_wt = min(p.grz_length * T(1e3) / dx, one(T))

    bmb_gl   = p.c_grz * bmb_ref - p.κ_grz * tf
    bmb_shelf = bmb_gl * f_shlf

    if mask == GroundingLine || mask == FloatingLine
        return (one(T) - grz_wt) * bmb_shelf + grz_wt * bmb_gl
    else
        return bmb_shelf
    end
end

# ============================================================
# In-place 2-D wrappers, dispatched on BMBMethod type
# ============================================================

"""
    compute_bmb!(bmb, now, par::BMBLinear, slope, dx)

Compute basal mass balance using the linear method and write into `bmb`.
If `par.γ_prime > 0` the exchange velocity is scaled locally by
`γ_prime · sin(θ)`.
"""
function compute_bmb!(
    bmb  ::AbstractMatrix{T},
    now  ::MarineShelfState{T},
    par  ::BMBLinear{T},
    slope::AbstractMatrix{T},
    dx   ::T,
) where T<:AbstractFloat
    ω_val = T(ω)
    use_slope = par.γ_prime > zero(T)
    @inbounds for j in axes(bmb, 2), i in axes(bmb, 1)
        γ = use_slope ? par.γ * par.γ_prime * slope[i,j] : par.γ
        bmb[i,j] = bmb_linear(now.tf_shlf[i,j], γ, ω_val)
    end
    return bmb
end

"""
    compute_bmb!(bmb, now, par::BMBQuad, slope, dx)

Compute basal mass balance using the local quadratic method.
"""
function compute_bmb!(
    bmb  ::AbstractMatrix{T},
    now  ::MarineShelfState{T},
    par  ::BMBQuad{T},
    slope::AbstractMatrix{T},
    dx   ::T,
) where T<:AbstractFloat
    ω_val = T(ω)
    use_slope = par.γ_prime > zero(T)
    @inbounds for j in axes(bmb, 2), i in axes(bmb, 1)
        γ = use_slope ? par.γ * par.γ_prime * slope[i,j] : par.γ
        bmb[i,j] = bmb_quad(now.tf_shlf[i,j], γ, ω_val)
    end
    return bmb
end

"""
    compute_bmb!(bmb, now, par::BMBQuadNL, slope, dx)

Compute basal mass balance using the non-local quadratic method.
Requires `now.tf_basin` to be pre-populated by `basin_mean!`.

When both `tf` and `tf_basin` are negative (cold conditions), the sign
of the result is flipped to allow ice accretion.
"""
function compute_bmb!(
    bmb  ::AbstractMatrix{T},
    now  ::MarineShelfState{T},
    par  ::BMBQuadNL{T},
    slope::AbstractMatrix{T},
    dx   ::T,
) where T<:AbstractFloat
    ω_val = T(ω)
    use_slope = par.γ_prime > zero(T)
    @inbounds for j in axes(bmb, 2), i in axes(bmb, 1)
        γ      = use_slope ? par.γ * par.γ_prime * slope[i,j] : par.γ
        tf     = now.tf_shlf[i,j]
        tf_bas = now.tf_basin[i,j]
        b      = bmb_quad_nl(tf, tf_bas, γ, ω_val)
        # Allow accretion when both local and basin TF are negative
        bmb[i,j] = (tf_bas ≤ zero(T) && tf ≤ zero(T)) ? -b : b
    end
    return bmb
end

"""
    compute_bmb!(bmb, now, par::BMBAnomaly, slope, dx)

Compute basal mass balance using the anomaly method.
`now.tf_shlf` must contain the temperature anomaly `dT_shlf`.
"""
function compute_bmb!(
    bmb  ::AbstractMatrix{T},
    now  ::MarineShelfState{T},
    par  ::BMBAnomaly{T},
    slope::AbstractMatrix{T},
    dx   ::T,
) where T<:AbstractFloat
    @inbounds for j in axes(bmb, 2), i in axes(bmb, 1)
        bmb[i,j] = bmb_anom(
            now.tf_shlf[i,j],
            now.bmb_ref[i,j],
            now.mask_ocn[i,j],
            par, dx,
        )
    end
    return bmb
end

"""
    compute_bmb!(bmb, now, par::BMBPico, slope, dx)

PICO placeholder — BMB is populated externally by the PICO sub-model and
copied into `now.bmb_shlf` before limits are applied.  This method is a
no-op; it exists only so dispatch succeeds when `par` is `BMBPico`.
"""
function compute_bmb!(
    bmb  ::AbstractMatrix{T},
    now  ::MarineShelfState{T},
    par  ::BMBPico{T},
    slope::AbstractMatrix{T},
    dx   ::T,
) where T<:AbstractFloat
    # BMB already populated externally (e.g. by pico_update! → now.bmb_shlf)
    return bmb
end

# ============================================================
# Post-processing limits
# ============================================================

"""
    apply_bmb_limits!(bmb, H_ice, mask_ocn, z_bed, z_sl, par)

Apply physical limits to `bmb` in-place after parameterisation:

1. Cap refreezing at `par.bmb_max`.
2. No accretion where there is no ice (`H_ice == 0`).
3. Apply `par.c_deep` melt rate at deep-ocean points deeper than
   `par.depth_deep`.
"""
function apply_bmb_limits!(
    bmb     ::AbstractMatrix{T},
    H_ice   ::AbstractMatrix{T},
    mask_ocn::AbstractMatrix{OceanMask},
    z_bed   ::AbstractMatrix{T},
    z_sl    ::AbstractMatrix{T},
    par     ::MarineShelfParams{T},
) where T<:AbstractFloat
    @inbounds for j in axes(bmb, 2), i in axes(bmb, 1)
        # 1. Cap refreezing
        if bmb[i,j] > par.bmb_max
            bmb[i,j] = par.bmb_max
        end

        # 2. No accretion without ice
        if bmb[i,j] > zero(T) && H_ice[i,j] == zero(T)
            bmb[i,j] = zero(T)
        end

        # 3. Deep-ocean forcing
        if mask_ocn[i,j] == DeepOcean &&
           (z_sl[i,j] - z_bed[i,j]) >= par.depth_deep
            bmb[i,j] = par.c_deep
        end
    end
    return bmb
end
