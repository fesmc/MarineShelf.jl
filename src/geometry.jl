# ============================================================
#  geometry.jl
#  Ice-shelf base geometry: elevation and slope.
# ============================================================

"""
    shelf_base_elevation!(z_base, H_ice, f_grnd, z_bed, z_sl)

Compute the ice-shelf base elevation relative to sea level [m] and write
it into the pre-allocated matrix `z_base`.

- Floating points:  `z_base = z_sl - H_ice · (ρ_ice/ρ_sw)`
- Grounded points:  `z_base = z_bed`  (used for geometry bookkeeping only)
"""
function shelf_base_elevation!(
    z_base::AbstractMatrix{T},
    H_ice ::AbstractMatrix{T},
    f_grnd::AbstractMatrix{T},
    z_bed ::AbstractMatrix{T},
    z_sl  ::AbstractMatrix{T},
) where T<:AbstractFloat
    ratio = T(ρ_ratio)
    @inbounds for j in axes(z_base, 2), i in axes(z_base, 1)
        if f_grnd[i,j] < one(T)
            # Floating (or partially floating) — shelf base
            z_base[i,j] = z_sl[i,j] - H_ice[i,j] * ratio
        else
            # Grounded — store bed elevation for completeness
            z_base[i,j] = z_bed[i,j]
        end
    end
    return z_base
end

"""
    shelf_base_slope!(slope, z_base, dx)

Compute `slope = |sin θ|` of the ice-shelf base at each grid point using
a centred finite difference, and write the result into `slope`.

    sin θ ≈ |Δz| / hypot(|Δz|, 2·dx)

where `Δz` is the mean of the x- and y-directed central differences.
Boundary points use one-sided differences (first/last column/row repeated).
"""
function shelf_base_slope!(
    slope ::AbstractMatrix{T},
    z_base::AbstractMatrix{T},
    dx    ::T,
) where T<:AbstractFloat
    nx, ny = size(slope)
    dx2 = 2 * dx
    @inbounds for j in 1:ny
        jm1 = max(j-1, 1); jp1 = min(j+1, ny)
        for i in 1:nx
            im1 = max(i-1, 1); ip1 = min(i+1, nx)
            dz = T(0.5) * (z_base[ip1,j] - z_base[im1,j]) +
                 T(0.5) * (z_base[i,jp1] - z_base[i,jm1])
            slope[i,j] = abs(dz) / sqrt(dz^2 + dx2^2)
        end
    end
    return slope
end
