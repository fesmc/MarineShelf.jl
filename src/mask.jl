# ============================================================
#  mask.jl
#  Ocean mask construction and open-ocean flood-fill.
# ============================================================

"""
    ocean_mask!(mask_ocn, f_grnd, mask_ocn_ref; find_ocean=true)

Classify each grid point and write the result into `mask_ocn`.

Pass 1 — per-point classification from `f_grnd`:
- `f_grnd == 1` with a fully-floating neighbour  → `GroundingLine`
- `f_grnd == 0` with a fully-grounded neighbour  → `FloatingLine`
- `f_grnd == 1` with no floating neighbour       → `Land`
- `f_grnd == 0` and `mask_ocn_ref == DeepOcean`  → `DeepOcean`
- otherwise                                       → `Ocean`

Pass 2 (optional, `find_ocean=true`) — flood-fill to mark points not
connected to the open ocean as `Lake`.
"""
function ocean_mask!(
    mask_ocn    ::AbstractMatrix{OceanMask},
    f_grnd      ::AbstractMatrix{T},
    mask_ocn_ref::AbstractMatrix{OceanMask};
    find_ocean  ::Bool = true,
) where T<:AbstractFloat
    nx, ny = size(mask_ocn)

    @inbounds for j in 1:ny
        jm1 = max(j-1,1); jp1 = min(j+1,ny)
        for i in 1:nx
            im1 = max(i-1,1); ip1 = min(i+1,nx)

            fij  = f_grnd[i,j]
            fn   = (f_grnd[im1,j], f_grnd[ip1,j],
                    f_grnd[i,jm1], f_grnd[i,jp1])
            any_float   = any(f -> f == zero(T), fn)
            any_grounded = any(f -> f == one(T), fn)

            if fij > zero(T) && any_float
                mask_ocn[i,j] = GroundingLine
            elseif fij == zero(T) && any_grounded
                mask_ocn[i,j] = FloatingLine
            elseif fij == one(T)
                mask_ocn[i,j] = Land
            elseif fij == zero(T) && mask_ocn_ref[i,j] == DeepOcean
                mask_ocn[i,j] = DeepOcean
            else
                mask_ocn[i,j] = Ocean
            end
        end
    end

    find_ocean && _find_open_ocean_and_lakes!(mask_ocn, mask_ocn_ref)
    return mask_ocn
end

# ------------------------------------------------------------
# Internal: flood-fill to separate open ocean from lakes
# ------------------------------------------------------------

"""
    _find_open_ocean_and_lakes!(mask_ocn, mask_ocn_ref)

Flood-fill from every point that is open ocean in `mask_ocn_ref`.
Any ocean or floating-line point not reachable from open ocean is
reclassified as `Lake`.  Grounding-line points with no open-ocean
neighbours are also reclassified as `Lake`.

Uses a working integer array to avoid allocating mask objects:
- `-2` = grounding line (deferred)
- `-1` = unknown / closed ocean
-  `0` = land
-  `1` = confirmed open ocean
"""
function _find_open_ocean_and_lakes!(
    mask_ocn    ::AbstractMatrix{OceanMask},
    mask_ocn_ref::AbstractMatrix{OceanMask},
)
    nx, ny = size(mask_ocn)
    work = zeros(Int8, nx, ny)

    @inbounds for j in 1:ny, i in 1:nx
        m = mask_ocn[i,j]
        if m == Land
            work[i,j] = Int8(0)
        elseif m == GroundingLine
            work[i,j] = Int8(-2)
        else
            # Ocean, DeepOcean, FloatingLine — initially unknown
            work[i,j] = Int8(-1)
        end
    end

    # Seed known open-ocean points from the reference mask
    @inbounds for j in 1:ny, i in 1:nx
        if work[i,j] == Int8(-1) &&
           (mask_ocn_ref[i,j] == Ocean || mask_ocn_ref[i,j] == DeepOcean)
            work[i,j] = Int8(1)
        end
    end

    # Iterative flood-fill (BFS-like single-pass repeated until convergence)
    max_iter = 1000
    for _ in 1:max_iter
        n_filled = 0
        @inbounds for j in 1:ny
            jm1 = max(j-1,1); jp1 = min(j+1,ny)
            for i in 1:nx
                im1 = max(i-1,1); ip1 = min(i+1,nx)
                if work[i,j] == Int8(1)
                    for (ni, nj) in ((im1,j),(ip1,j),(i,jm1),(i,jp1))
                        if work[ni,nj] == Int8(-1)
                            work[ni,nj] = Int8(1)
                            n_filled += 1
                        end
                    end
                end
            end
        end
        n_filled == 0 && break
    end

    # Reclassify unreached ocean/floating-line points as lakes
    @inbounds for j in 1:ny, i in 1:nx
        m = mask_ocn[i,j]
        if work[i,j] == Int8(-1) && (m == Ocean || m == FloatingLine)
            mask_ocn[i,j] = Lake
        end
    end

    # Grounding-line points with no open-ocean neighbour → lake
    @inbounds for j in 1:ny
        jm1 = max(j-1,1); jp1 = min(j+1,ny)
        for i in 1:nx
            im1 = max(i-1,1); ip1 = min(i+1,nx)
            if work[i,j] == Int8(-2)
                has_open = work[im1,j] == Int8(1) || work[ip1,j] == Int8(1) ||
                           work[i,jm1] == Int8(1) || work[i,jp1] == Int8(1)
                has_open || (mask_ocn[i,j] = Lake)
            end
        end
    end

    return mask_ocn
end

"""
    build_mask_ocn_ref!(mask_ocn_ref, regions; domain=nothing)

Populate the static reference mask from a `regions` field.

Rules common to all domains:
- `regions == 1` or `regions == 2` → `Ocean`
- everything else                  → `Land`

Domain-specific rules (pass `domain` keyword):
- `"Greenland"`:  all non-1.3 regions → `DeepOcean`
- `"Antarctica"`: all non-2.11 regions → `DeepOcean`
- other / nothing: deep ocean everywhere outside land

These match the legacy Fortran behaviour; override by mutating
`mask_ocn_ref` directly after calling this function.
"""
function build_mask_ocn_ref!(
    mask_ocn_ref::AbstractMatrix{OceanMask},
    regions     ::AbstractMatrix{T};
    domain      ::Union{Nothing,AbstractString} = nothing,
) where T<:AbstractFloat
    fill!(mask_ocn_ref, Land)
    @inbounds for j in axes(regions, 2), i in axes(regions, 1)
        r = regions[i,j]
        if r == one(T) || r == 2*one(T)
            mask_ocn_ref[i,j] = Ocean
        end
    end

    if domain == "Greenland"
        @inbounds for j in axes(regions, 2), i in axes(regions, 1)
            if regions[i,j] != T(1.3)
                mask_ocn_ref[i,j] = DeepOcean
            end
        end
    elseif domain == "Antarctica"
        @inbounds for j in axes(regions, 2), i in axes(regions, 1)
            if regions[i,j] != T(2.11)
                mask_ocn_ref[i,j] = DeepOcean
            end
        end
    else
        # Generic: deep ocean everywhere (c_deep applied everywhere)
        @inbounds for j in axes(regions, 2), i in axes(regions, 1)
            if mask_ocn_ref[i,j] != Land
                mask_ocn_ref[i,j] = DeepOcean
            end
        end
    end

    return mask_ocn_ref
end
