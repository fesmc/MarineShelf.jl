# ============================================================
#  io.jl
#  NetCDF restart I/O using NCDatasets.jl.
#
#  Design:
#  - `write_restart` appends a time slice to an existing file or
#    creates the file on first call.
#  - `read_restart!` loads the first (or last) time slice back
#    into a MarineShelfState.
#  - The file stores only the fields needed to warm-restart a run;
#    parameters are not stored (load those from your namelist/TOML).
# ============================================================

using NCDatasets

# Fields written / read during restart (name → accessor)
const _RESTART_FIELDS = (
    # name            accessor in MarineShelfState
    ("bmb_shlf",       s -> s.bmb_shlf),
    ("bmb_ref",        s -> s.bmb_ref),
    ("bmb_corr",       s -> s.bmb_corr),
    ("T_shlf",         s -> s.T_shlf),
    ("S_shlf",         s -> s.S_shlf),
    ("dT_shlf",        s -> s.dT_shlf),
    ("T_fp_shlf",      s -> s.T_fp_shlf),
    ("tf_shlf",        s -> s.tf_shlf),
    ("tf_basin",       s -> s.tf_basin),
    ("tf_corr",        s -> s.tf_corr),
    ("tf_corr_basin",  s -> s.tf_corr_basin),
    ("T_shlf_basin",   s -> s.T_shlf_basin),
    ("S_shlf_basin",   s -> s.S_shlf_basin),
    ("z_base",         s -> s.z_base),
    ("slope_base",     s -> s.slope_base),
)

"""
    write_restart(path, mshlf, xc, yc, time; append=true)

Write all restart fields from `mshlf.now` to a NetCDF file at `path`.

- On first call (or `append=false`) the file is created and dimensions
  are defined.
- On subsequent calls the time dimension is extended.

`xc`, `yc` are 1-D coordinate vectors [m]; they are written as [km].
`time` is the model time [yr].
"""
function write_restart(
    path  ::AbstractString,
    mshlf ::MarineShelf{T},
    xc    ::AbstractVector{T},
    yc    ::AbstractVector{T},
    time  ::T;
    append::Bool = true,
) where T<:AbstractFloat
    now = mshlf.now
    nx, ny = size(now.bmb_shlf)

    create_new = !isfile(path) || !append

    NCDataset(path, create_new ? "c" : "a") do ds
        if create_new
            # Dimensions
            defDim(ds, "xc",   nx)
            defDim(ds, "yc",   ny)
            defDim(ds, "time", Inf)   # unlimited

            # Coordinate variables
            xv = defVar(ds, "xc",   Float64, ("xc",),   attrib=["units"=>"km"])
            yv = defVar(ds, "yc",   Float64, ("yc",),   attrib=["units"=>"km"])
            tv = defVar(ds, "time", Float64, ("time",),  attrib=["units"=>"years"])
            xv[:] = xc .* 1e-3
            yv[:] = yc .* 1e-3

            # Data variables
            for (name, _) in _RESTART_FIELDS
                # skip duplicates from the tuple list
                haskey(ds, name) && continue
                defVar(ds, name, Float32, ("xc","yc","time"))
            end
            # Masks stored as integers
            defVar(ds, "mask_ocn_ref", Int8, ("xc","yc","time"))
            defVar(ds, "mask_ocn",     Int8, ("xc","yc","time"))
        end

        # Time index to write
        nt = ds.dim["time"]
        n  = create_new ? 1 : nt + 1

        ds["time"][n] = Float64(time)

        for (name, acc) in _RESTART_FIELDS
            haskey(ds, name) || continue
            ds[name][:,:,n] = Float32.(acc(now))
        end

        ds["mask_ocn_ref"][:,:,n] = Int8.(Int.(now.mask_ocn_ref))
        ds["mask_ocn"][:,:,n]     = Int8.(Int.(now.mask_ocn))
    end

    return nothing
end

"""
    read_restart!(mshlf, path; time_index=1)

Read restart fields from `path` into `mshlf.now`.
`time_index` selects which time slice to load (1-based; use `nothing`
to load the last available slice).
"""
function read_restart!(
    mshlf     ::MarineShelf{T},
    path      ::AbstractString;
    time_index::Union{Int,Nothing} = 1,
) where T<:AbstractFloat
    now = mshlf.now

    NCDataset(path, "r") do ds
        n = time_index === nothing ? ds.dim["time"] : time_index

        for (name, acc) in _RESTART_FIELDS
            haskey(ds, name) || continue
            field = acc(now)
            field .= T.(ds[name][:,:,n])
        end

        raw_ref = ds["mask_ocn_ref"][:,:,n]
        raw_ocn = ds["mask_ocn"][:,:,n]
        @inbounds for j in axes(now.mask_ocn_ref, 2),
                      i in axes(now.mask_ocn_ref, 1)
            now.mask_ocn_ref[i,j] = OceanMask(raw_ref[i,j])
            now.mask_ocn[i,j]     = OceanMask(raw_ocn[i,j])
        end
    end

    return mshlf
end
