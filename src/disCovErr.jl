module disCovErr

    export cov_construct
    export per_star
    export per_star_stats
    export prelim_infill!
    export add_sky_noise!
    export gen_mask_staticPSF!

    include("utils.jl")
    include("decam.jl")
end
