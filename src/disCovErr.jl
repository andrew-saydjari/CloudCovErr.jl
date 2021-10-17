module disCovErr

    export cov_construct
    export per_star
    export per_star_stats
    export prelim_infill!
    export add_sky_noise!
    export gen_mask_staticPSF!
    export condCovEst_wdiag
    export stamp_cutter
    export read_decam
    export read_crowdsource
    export proc_ccd
    export gen_mask_staticPSF!
    export prelim_infill!
    export add_sky_noise!
    export load_psfmodel_cs
    export stamp_cutter
    export gen_pix_mask
    export condCovEst_wdiag

    include("utils.jl")
    include("decam.jl")
    include("decaps2.jl")
end
