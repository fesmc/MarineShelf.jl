# ============================================================
#  extrapolation.jl
#  Horizontal extrapolation of 3-D ocean fields into
#  ice-shelf cavities.
# ============================================================

"""
    ocn_variable_extrapolation!(var, H_ice, basin_mask, depth, z_bed)

Extrapolate a 3-D ocean field `var[nx, ny, nz]` into ice-shelf cavities
and other regions with no direct ocean data (identified by missing values
or by the presence of overlying ice).

The algorithm (matching `ocn_variable_extrapolation` in the Fortran module):

1. Mark all points with ice (`H_ice > 0`) or missing values as "cavity"
   points.  Identify the ocean border — the open-ocean points immediately
   adjacent to cavity points.
2. For each basin and depth level, average `var` over the border points
   whose bed depth exceeds the current level depth.
3. Fill every cavity point with the basin-and-depth-appropriate mean.
   If a depth level has no border points, carry the mean from the level
   above (vertical extrapolation to handle sills).

## Arguments
- `var`        — 3-D ocean field `[nx, ny, nz]`, modified in-place.
- `H_ice`      — ice thickness `[nx, ny]` [m]; points with `H_ice > 0`
                  are treated as cavities.
- `basin_mask` — integer basin index `[nx, ny]`; basin `0` is ignored
                  (outside the continental-shelf break).
- `depth`      — 1-D depth coordinate `[nz]` [m], positive downward.
- `z_bed`      — bed elevation `[nx, ny]` [m], positive upward (negative
                  below sea level).

## Missing-value convention
The missing-value sentinel is `mv = -9999.0`.  Points where
`var[:,:,1] == mv` are also treated as cavities.
"""
function ocn_variable_extrapolation!(
    var       ::AbstractArray{T,3},
    H_ice     ::AbstractMatrix{T},
    basin_mask::AbstractMatrix,
    depth     ::AbstractVector{T},
    z_bed     ::AbstractMatrix{T},
) where T<:AbstractFloat

    nx, ny, nz = size(var)
    nb = Int(maximum(basin_mask))

    mv = T(-9999)

    # ---- Step 1: classify points as cavity (2), border (1), or ocean (0)
    mask_border = zeros(T, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        if H_ice[i,j] > zero(T) || var[i,j,1] == mv
            mask_border[i,j] = T(2)
        end
    end

    @inbounds for j in 2:ny-1, i in 2:nx-1
        if mask_border[i,j] == zero(T)
            if mask_border[i+1,j] == T(2) || mask_border[i-1,j] == T(2) ||
               mask_border[i,j+1] == T(2) || mask_border[i,j-1] == T(2)
                mask_border[i,j] = one(T)
            end
        end
    end

    # ---- Step 2: compute basin-and-depth-mean over border points
    mean_var_basin = zeros(T, nb, nz)
    npts_basin     = zeros(T, nb, nz)

    @inbounds for k in 1:nz
        for j in 1:ny, i in 1:nx
            if mask_border[i,j] == one(T) && z_bed[i,j] < depth[k]
                b = Int(basin_mask[i,j])
                b > 0 || continue
                mean_var_basin[b,k] += var[i,j,k]
                npts_basin[b,k]     += one(T)
            end
        end
    end

    @. mean_var_basin = mean_var_basin / (npts_basin + T(1e-8))

    # Vertical extrapolation: if no border points at a level, use level above
    @inbounds for b in 1:nb, k in 2:nz
        if npts_basin[b,k] == zero(T)
            mean_var_basin[b,k] = mean_var_basin[b,k-1]
        end
    end

    # ---- Step 3: fill cavity points
    @inbounds for k in 1:nz
        for j in 1:ny, i in 1:nx
            if mask_border[i,j] == T(2)
                b = Int(basin_mask[i,j])
                b > 0 || continue
                var[i,j,k] = mean_var_basin[b,k]
            end
        end
    end

    return var
end
