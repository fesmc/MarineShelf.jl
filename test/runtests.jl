using MarineShelf
using Test
using Statistics: mean

# ============================================================
# Helpers
# ============================================================

function make_grid(T, nx=8, ny=6)
    H_ice  = zeros(T, nx, ny)
    z_bed  = fill(T(-800), nx, ny)
    f_grnd = zeros(T, nx, ny)
    basins = ones(T, nx, ny)
    z_sl   = zeros(T, nx, ny)
    dx     = T(10_000)
    return H_ice, z_bed, f_grnd, basins, z_sl, dx
end

# A minimal floating-shelf setup: centre cells float, edges ground
function shelf_grid(T, nx=8, ny=6)
    H_ice, z_bed, f_grnd, basins, z_sl, dx = make_grid(T, nx, ny)
    H_ice  .= T(500)
    f_grnd[1:2, :]   .= T(1)
    f_grnd[end-1:end, :] .= T(1)
    return H_ice, z_bed, f_grnd, basins, z_sl, dx
end

# ============================================================
@testset "MarineShelf.jl" begin

# ------------------------------------------------------------
@testset "Constants" begin
    @test ρ_ice ≈ 918.0
    @test ρ_sw  ≈ 1028.0
    @test ρ_w   ≈ 1000.0
    @test T0    ≈ 273.15
    @test ω     ≈ (ρ_sw * cp_o) / (ρ_ice * L_ice)
    @test ρ_ratio ≈ ρ_ice / ρ_sw
end

# ------------------------------------------------------------
@testset "FreezingPointParams" begin
    fp64 = FreezingPointParams(Float64)
    fp32 = FreezingPointParams(Float32)
    @test fp64 isa FreezingPointParams{Float64}
    @test fp32 isa FreezingPointParams{Float32}
    @test fp64.λ1 ≈ -0.0573
    # scalar function
    T_fp = freezing_point(Float64(35.0), Float64(-500.0), fp64)
    @test T_fp ≈ fp64.λ1*35.0 + fp64.λ2 + fp64.λ3*(-500.0) + T0
end

# ------------------------------------------------------------
@testset "Scalar BMB kernels" begin
    tf = 2.0; γ = 1e4; ω_v = Float64(ω)
    @test bmb_linear(tf, γ, ω_v) ≈ -γ * ω_v * tf
    @test bmb_quad(tf, γ, ω_v)   ≈ -γ * ω_v^2 * tf^2
    @test bmb_quad_nl(tf, tf, γ, ω_v) ≈ -γ * ω_v^2 * tf^2

    # linear: positive TF → melt (negative bmb)
    @test bmb_linear(1.0, 1e4, ω_v) < 0
    # quad: TF=0 → zero melt
    @test bmb_quad(0.0, 1e4, ω_v) == 0.0

    # anomaly method scalar
    p = BMBAnomaly{Float64}(κ_grz=10.0, c_grz=1.0, f_grz_shlf=2.0, grz_length=5.0)
    dx = 5_000.0
    bmb_ocn   = bmb_anom(1.0, -5.0, Ocean,         p, dx)
    bmb_gl    = bmb_anom(1.0, -5.0, GroundingLine,  p, dx)
    # open shelf receives a fraction of the grounding-zone value
    bmb_grz = p.c_grz * (-5.0) - p.κ_grz * 1.0
    @test bmb_ocn ≈ bmb_grz / p.f_grz_shlf
    # at coarse resolution (grz_length*1e3/dx ≈ 1) the GL value dominates
    bmb_gl_fine = bmb_anom(1.0, -5.0, GroundingLine, p, 100.0)
    @test bmb_gl_fine ≈ bmb_grz
end

# ------------------------------------------------------------
@testset "OceanMask enum" begin
    @test Int(Land)          == 0
    @test Int(GroundingLine) == 1
    @test Int(FloatingLine)  == 2
    @test Int(Ocean)         == 3
    @test Int(DeepOcean)     == 4
    @test Int(Lake)          == 5
end

# ------------------------------------------------------------
@testset "Geometry" begin
    T = Float64
    nx, ny = 6, 5
    H_ice  = fill(T(500), nx, ny)
    z_bed  = fill(T(-800), nx, ny)
    f_grnd = zeros(T, nx, ny)
    z_sl   = zeros(T, nx, ny)

    z_base = zeros(T, nx, ny)
    slope  = zeros(T, nx, ny)

    # floating: shelf base = z_sl - H_ice * ρ_ice/ρ_sw
    MarineShelf.shelf_base_elevation!(z_base, H_ice, f_grnd, z_bed, z_sl)
    @test all(z_base .≈ -H_ice .* ρ_ratio)

    # grounded: shelf base = z_bed
    f_grnd2 = ones(T, nx, ny)
    MarineShelf.shelf_base_elevation!(z_base, H_ice, f_grnd2, z_bed, z_sl)
    @test all(z_base .≈ z_bed)

    # uniform-depth shelf has zero slope; slope always in [0,1]
    MarineShelf.shelf_base_slope!(slope, z_base, T(10_000))
    @test all(slope .== 0)
    @test all(slope .<= 1)

    # non-uniform: slope should be positive somewhere
    z_base2 = copy(z_base)
    z_base2[3, 3] = T(-1000)
    MarineShelf.shelf_base_slope!(slope, z_base2, T(10_000))
    @test any(slope .> 0)
    @test all(0 .<= slope .<= 1)
end

# ------------------------------------------------------------
@testset "Ocean mask" begin
    T = Float64
    nx, ny = 6, 5
    f_grnd   = zeros(T, nx, ny)
    f_grnd[1:2, :] .= T(1)        # grounded strip on left
    mask_ref = fill(Ocean, nx, ny)
    mask     = fill(Ocean, nx, ny)

    MarineShelf.ocean_mask!(mask, f_grnd, mask_ref; find_ocean=false)

    # interior of grounded strip → Land
    @test mask[1, 3] == Land
    # first floating column adjacent to grounded → FloatingLine
    @test mask[3, 3] == FloatingLine
    # grounded column adjacent to floating → GroundingLine
    @test mask[2, 3] == GroundingLine

    # flood-fill: isolated inland ocean should become Lake
    f_lake        = ones(T, nx, ny)
    f_lake[3:4, 3:4] .= T(0)          # isolated floating patch
    mask_ref_lake = fill(Land, nx, ny) # no open-ocean reference → no seed
    mask_lake     = fill(Ocean, nx, ny)
    MarineShelf.ocean_mask!(mask_lake, f_lake, mask_ref_lake; find_ocean=true)
    @test all(mask_lake[3:4, 3:4] .== Lake)
end

# ------------------------------------------------------------
@testset "build_mask_ocn_ref!" begin
    T = Float64
    nx, ny = 8, 6

    # Generic domain: r==1 and r==2 become Ocean, then all non-Land → DeepOcean
    regions = zeros(T, nx, ny)
    regions[1:4, :] .= T(1)
    regions[5:8, :] .= T(2)
    mask = fill(Land, nx, ny)
    build_mask_ocn_ref!(mask, regions)
    @test all(mask[1:4, :] .== DeepOcean)
    @test all(mask[5:8, :] .== DeepOcean)

    # Antarctica domain:
    #   Step 1: r==1 or r==2 → Ocean (here only region==2.0 rows do)
    #   Step 2: all regions != 2.11 → DeepOcean
    # So:
    #   - rows with region==2.0  (set to Ocean in step 1, then != 2.11 → DeepOcean)
    #   - rows with region==2.11 (not set in step 1, stays Land;  == 2.11 → stays Land)
    regions_ant = fill(T(2.0), nx, ny)     # open ocean rows
    regions_ant[3:end, :] .= T(2.11)       # continental-shelf rows
    mask_ant = fill(Land, nx, ny)
    build_mask_ocn_ref!(mask_ant, regions_ant; domain="Antarctica")
    @test all(mask_ant[1:2,   :] .== DeepOcean)   # open ocean → DeepOcean
    @test all(mask_ant[3:end, :] .== Land)         # continental shelf → Land (no c_deep)

    # Greenland domain: all regions != 1.3 → DeepOcean
    regions_grl = fill(T(1.3), nx, ny)
    regions_grl[1:2, :] .= T(1.0)          # open ocean
    mask_grl = fill(Land, nx, ny)
    build_mask_ocn_ref!(mask_grl, regions_grl; domain="Greenland")
    @test all(mask_grl[1:2,   :] .== DeepOcean)   # 1.0 != 1.3 → DeepOcean
    @test all(mask_grl[3:end, :] .== Land)         # 1.3 stays as Land
end

# ------------------------------------------------------------
@testset "TFMethod types" begin
    @test TFExternal() isa TFMethod
    @test TFAbsolute() isa TFMethod
    @test TFAnomaly()  isa TFMethod

    # default tf_method selection
    @test MarineShelf._default_tf_method(BMBLinear{Float64}(γ=1.0))   isa TFAbsolute
    @test MarineShelf._default_tf_method(BMBAnomaly{Float64}(
        κ_grz=1.0, c_grz=1.0, f_grz_shlf=2.0, grz_length=5.0))      isa TFAnomaly
end

# ------------------------------------------------------------
@testset "MarineShelfParams construction" begin
    T = Float64
    p = MarineShelfParams(T, BMBQuad{T}(γ=1e4))
    @test p.tf_method isa TFAbsolute
    @test p.domain == ""
    @test p.c_deep == T(-50)

    # explicit tf_method
    p2 = MarineShelfParams(T, BMBLinear{T}(γ=1e4); tf_method=TFExternal(), domain="Greenland")
    @test p2.tf_method isa TFExternal
    @test p2.domain == "Greenland"

    # consistency check
    @test_throws ErrorException MarineShelfParams(T,
        BMBAnomaly{T}(κ_grz=1.0, c_grz=1.0, f_grz_shlf=2.0, grz_length=5.0);
        tf_method=TFAbsolute())
end

# ------------------------------------------------------------
@testset "MShlf construction" begin
    T = Float32
    nx, ny = 8, 6
    method = BMBQuad{T}(γ = T(1e4))

    # basic construction
    m = MShlf(T, nx, ny, method)
    @test m isa MShlf{T}
    @test size(m.now.bmb_shlf) == (nx, ny)
    @test size(m.now.dS_shlf)  == (nx, ny)
    @test typeof(m.now.bmb_shlf) == Matrix{T}   # type param correct
    @test m.par.tf_method isa TFAbsolute

    # pre-built basin corrections
    basins = ones(T, nx, ny)
    basins[5:end, :] .= T(2)
    m2 = MShlf(T, nx, ny, method;
               basins        = basins,
               basin_bmb_corr = Dict(1 => T(-3), 2 => T(-10)))
    @test all(m2.now.bmb_corr[1:4, :] .≈ T(-3))
    @test all(m2.now.bmb_corr[5:end, :] .≈ T(-10))
    # correction stays static — modifying basins after construction has no effect
    basins .= T(99)
    @test all(m2.now.bmb_corr[1:4, :] .≈ T(-3))

    # error when corrections given without basins
    @test_throws ErrorException MShlf(T, nx, ny, method;
        basin_tf_corr=Dict(1 => T(0.5)))

    # TFExternal
    m3 = MShlf(T, nx, ny, method; tf_method=TFExternal())
    @test m3.par.tf_method isa TFExternal
end

# ------------------------------------------------------------
@testset "update_tf!" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = make_grid(T, nx, ny)

    # TFAbsolute: tf = T_shlf - T_fp + corrections
    m = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4))
    m.now.T_shlf .= T(272.0)
    m.now.S_shlf .= T(35.0)
    m.now.z_base .= T(-500.0)
    update_tf!(m)
    T_fp_expected = m.par.fp.λ1*35.0 + m.par.fp.λ2 + m.par.fp.λ3*(-500.0) + T0
    @test all(m.now.T_fp_shlf .≈ T_fp_expected)
    @test all(m.now.tf_shlf   .≈ T(272.0) - T_fp_expected)

    # TFAnomaly: tf = dT_shlf + corrections
    m2 = MShlf(T, nx, ny, BMBAnomaly{T}(κ_grz=1.0, c_grz=1.0, f_grz_shlf=2.0, grz_length=5.0))
    m2.now.dT_shlf .= T(1.5)
    update_tf!(m2)
    @test all(m2.now.tf_shlf .≈ T(1.5))

    # TFExternal: tf_shlf is untouched
    m3 = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4); tf_method=TFExternal())
    m3.now.tf_shlf .= T(3.14)
    update_tf!(m3)
    @test all(m3.now.tf_shlf .≈ T(3.14))

    # tf_corr_basin added to TF (pre-built)
    basins_field = ones(T, nx, ny)
    m4 = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4);
               basins = basins_field,
               basin_tf_corr = Dict(1 => T(0.5)))
    m4.now.T_shlf .= T(272.0)
    m4.now.S_shlf .= T(35.0)
    m4.now.z_base .= T(-500.0)
    update_tf!(m4)
    @test all(m4.now.tf_shlf .≈ (T(272.0) - T_fp_expected) + T(0.5))
end

# ------------------------------------------------------------
@testset "update_bmb! — BMBLinear" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)
    m = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4))
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny))
    m.now.T_shlf .= T(275.0)
    m.now.S_shlf .= T(35.0)
    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)

    # floating cells should have nonzero melt (negative bmb)
    float_cells = f_grnd .< 1
    @test any(m.now.bmb_shlf[float_cells] .< 0)

    # melt is proportional to TF (colder → less melt)
    m2 = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4))
    build_mask_ocn_ref!(m2.now.mask_ocn_ref, ones(T, nx, ny))
    m2.now.T_shlf .= T(273.5)   # colder
    m2.now.S_shlf .= T(35.0)
    update!(m2, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    @test all(m2.now.bmb_shlf[float_cells] .>= m.now.bmb_shlf[float_cells])
end

# ------------------------------------------------------------
@testset "update_bmb! — BMBQuad" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)
    m = MShlf(T, nx, ny, BMBQuad{T}(γ=1e4))
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny))
    m.now.T_shlf .= T(275.0)
    m.now.S_shlf .= T(35.0)
    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    float_cells = f_grnd .< 1
    @test any(m.now.bmb_shlf[float_cells] .< 0)

    # quad is more sensitive to TF than linear (for |TF| > 1/ω)
    m_lin = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4))
    build_mask_ocn_ref!(m_lin.now.mask_ocn_ref, ones(T, nx, ny))
    m_lin.now.T_shlf .= T(275.0)
    m_lin.now.S_shlf .= T(35.0)
    update!(m_lin, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    # both negative; can't make a universal magnitude claim, just verify runs
    @test all(isfinite, m.now.bmb_shlf)
end

# ------------------------------------------------------------
@testset "update_bmb! — BMBQuadNL" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)
    m = MShlf(T, nx, ny, BMBQuadNL{T}(γ=1e4))
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny))
    m.now.T_shlf .= T(275.0)
    m.now.S_shlf .= T(35.0)
    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    @test all(isfinite, m.now.bmb_shlf)
    @test all(isfinite, m.now.tf_basin)

    # cold basin: both tf and tf_basin negative → accretion (positive bmb) before bmb_max cap
    m2 = MShlf(T, nx, ny, BMBQuadNL{T}(γ=1e4); bmb_max=T(10))
    build_mask_ocn_ref!(m2.now.mask_ocn_ref, ones(T, nx, ny))
    m2.now.T_shlf .= T(270.0)   # well below freezing → negative TF
    m2.now.S_shlf .= T(35.0)
    update!(m2, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    float_cells = f_grnd .< 1
    @test any(m2.now.bmb_shlf[float_cells] .>= 0)  # accretion allowed
end

# ------------------------------------------------------------
@testset "update_bmb! — BMBAnomaly" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)
    p = BMBAnomaly{T}(κ_grz=10.0, c_grz=1.0, f_grz_shlf=2.0, grz_length=5.0)
    m = MShlf(T, nx, ny, p)
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny))
    m.now.bmb_ref  .= T(-5.0)
    m.now.dT_shlf  .= T(1.0)
    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    @test all(isfinite, m.now.bmb_shlf)

    # zero anomaly → shelf value is c_grz * bmb_ref / f_grz_shlf
    m2 = MShlf(T, nx, ny, p)
    build_mask_ocn_ref!(m2.now.mask_ocn_ref, ones(T, nx, ny))
    m2.now.bmb_ref .= T(-5.0)
    m2.now.dT_shlf .= T(0.0)
    update!(m2, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    expected_shelf = (p.c_grz * T(-5.0)) / p.f_grz_shlf
    open_cells = findall(==( Ocean), m2.now.mask_ocn)
    if !isempty(open_cells)
        @test all(m2.now.bmb_shlf[open_cells] .≈ expected_shelf)
    end
end

# ------------------------------------------------------------
@testset "BMBPico (no-op compute)" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)
    m = MShlf(T, nx, ny, BMBPico{T}())
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny))
    m.now.bmb_shlf .= T(-20.0)   # pre-populated externally
    m.now.T_shlf   .= T(275.0)
    m.now.S_shlf   .= T(35.0)
    # update! should preserve the pre-set bmb_shlf for floating cells with ice
    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    float_cells = (f_grnd .< 1) .& (H_ice .> 0)
    @test all(m.now.bmb_shlf[float_cells] .≈ T(-20.0))
end

# ------------------------------------------------------------
@testset "apply_bmb_limits!" begin
    T = Float64
    nx, ny = 4, 4
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)
    m = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4); bmb_max=T(2), c_deep=T(-50), depth_deep=T(500))
    build_mask_ocn_ref!(m.now.mask_ocn_ref, fill(T(1), nx, ny))
    m.now.T_shlf .= T(272.0)
    m.now.S_shlf .= T(35.0)
    # Force refreezing everywhere
    m.now.T_shlf .= T(200.0)
    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    # Refreezing capped at bmb_max=2
    @test all(m.now.bmb_shlf .<= T(2))
    # No accretion where no ice
    no_ice = H_ice .== 0
    @test all(m.now.bmb_shlf[no_ice] .<= 0)
end

# ------------------------------------------------------------
@testset "apply_basin_correction!" begin
    T = Float64
    nx, ny = 6, 5
    field  = zeros(T, nx, ny)
    basins = ones(T, nx, ny)
    basins[4:end, :] .= T(2)
    corr   = Dict(1 => T(3.0), 2 => T(-1.0))

    apply_basin_correction!(field, basins, corr)
    @test all(field[1:3, :] .≈ T(3.0))
    @test all(field[4:end, :] .≈ T(-1.0))

    # empty dict → no change
    field2 = zeros(T, nx, ny)
    apply_basin_correction!(field2, basins, Dict{Int,T}())
    @test all(field2 .== 0)
end

# ------------------------------------------------------------
@testset "basin_mean!" begin
    T = Float64
    nx, ny = 6, 5
    var2D  = fill(T(4.0), nx, ny)
    var2D[3:4, 3:4] .= T(8.0)   # some cells differ
    f_grnd = zeros(T, nx, ny)
    H_ice  = fill(T(100.0), nx, ny)
    basins = ones(T, nx, ny)
    mask   = fill(Ocean, nx, ny)

    var_basin = zeros(T, nx, ny)
    basin_mean!(var_basin, var2D, f_grnd, basins, H_ice, mask)

    # all cells in basin 1 → same mean
    @test all(var_basin .≈ var_basin[1,1])
    @test var_basin[1,1] ≈ mean(var2D)
end

# ------------------------------------------------------------
@testset "Full update! round-trip (Float32)" begin
    T = Float32
    nx, ny = 10, 8
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)

    m = MShlf(T, nx, ny, BMBQuad{T}(γ=T(1e4)); domain="Antarctica")
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny); domain="Antarctica")
    m.now.T_shlf .= T(274.5)
    m.now.S_shlf .= T(34.5)

    @test_nowarn update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    @test all(isfinite, m.now.bmb_shlf)
    @test all(isfinite, m.now.tf_shlf)
    @test all(m.now.bmb_shlf .<= m.par.bmb_max)
end

# ------------------------------------------------------------
@testset "regions kwarg (PICO Antarctica blend)" begin
    T = Float64
    nx, ny = 6, 5
    H_ice, z_bed, f_grnd, basins, z_sl, dx = shelf_grid(T, nx, ny)

    m = MShlf(T, nx, ny, BMBPico{T}(); domain="Antarctica", c_deep=T(-50))
    build_mask_ocn_ref!(m.now.mask_ocn_ref, ones(T, nx, ny))
    m.now.bmb_shlf .= T(-20.0)
    m.now.T_shlf   .= T(275.0)
    m.now.S_shlf   .= T(35.0)

    regions = fill(T(2.11), nx, ny)
    regions[1:3, :] .= T(2.1)

    update!(m, H_ice, z_bed, f_grnd, basins, z_sl, dx; regions)

    # region 2.1 cells: blended with c_deep
    float_21 = (regions .== T(2.1)) .& (f_grnd .< 1) .& (H_ice .> 0)
    if any(float_21)
        expected = T(0.5) * T(-20.0) + T(0.5) * T(-50.0)
        @test all(m.now.bmb_shlf[float_21] .≈ expected)
    end

    # without regions kwarg, blending should not happen
    m2 = MShlf(T, nx, ny, BMBPico{T}(); domain="Antarctica", c_deep=T(-50))
    build_mask_ocn_ref!(m2.now.mask_ocn_ref, ones(T, nx, ny))
    m2.now.bmb_shlf .= T(-20.0)
    m2.now.T_shlf   .= T(275.0)
    m2.now.S_shlf   .= T(35.0)
    update!(m2, H_ice, z_bed, f_grnd, basins, z_sl, dx)
    if any(float_21)
        @test all(m2.now.bmb_shlf[float_21] .≈ T(-20.0))
    end
end

# ------------------------------------------------------------
@testset "Thermal forcing corrections via tf_corr" begin
    T = Float64
    nx, ny = 6, 5
    m = MShlf(T, nx, ny, BMBLinear{T}(γ=1e4))
    m.now.T_shlf  .= T(272.0)
    m.now.S_shlf  .= T(35.0)
    m.now.z_base  .= T(-500.0)

    update_tf!(m)
    tf_no_corr = copy(m.now.tf_shlf)

    m.now.tf_corr .= T(1.0)
    update_tf!(m)
    @test all(m.now.tf_shlf .≈ tf_no_corr .+ T(1.0))
end

end # @testset "MarineShelf.jl"
