using disCovErr
using ArgParse
using FITSIO #this should not be necessary

function parse_commandline()
    s=ArgParseSettings()
    @add_arg_table s begin
        "base"
            required=true
            arg_type=String
            help="base directory containing raw exposure image files"
        "date"
            required=true
            arg_type=String
            help="file name date (exposure datetime)"
        "filt"
            required=true
            arg_type=String
            help="file name filter"
        "vers"
            required=true
            arg_type=String
            help="file name version"
        "basecat"
            required=true
            arg_type=String
            help="base directory containing crowdsource cat files"
        "--thr", "-t"
            help = "threshold for psf-based masking of the residuals (int for ease)"
            arg_type = Int
            default = 20
        "--Np", "-p"
            help = "sidelength size of spatial covariance matrix (must be int, odd)"
            arg_type = Int
            default = 33
    end
end

function run_wrapper()
    parg = parse_commandline()
    proc_ccd(parg["base"],parg["date"],parg["filt"],parg["vers"],
        parg["basecat"],"N14",thr=parg["thr"],Np=parg["Np"])
end

# should we be prellocating outside this subfunction?
function proc_ccd(base,date,filt,vers,basecat,ccd;thr=20,Np=33)
    # loads from disk
    ref_im, w_im, d_im = read_decam(base,date,filt,vers,ccd)
    (sx, sy) = size(ref_im)
    x_stars, y_stars, flux_stars, decapsid, gain, mod_im, sky_im, w = read_crowdsource(basecat,date,filt,vers,ccd)

    psfmodel = load_psfmodel_cs(basecat,date,filt,vers,ccd)
    psfstatic = psfmodel(sx÷2,sy÷2,511)

    # mask bad camera pixels/cosmic rays, then mask out star centers
    bmaskd = (d_im .!= 0)
    gen_mask_staticPSF!(bmaskd, psfstatic, x_stars, y_stars, flux_stars; thr=thr)

    testim = copy(mod_im .- ref_im)
    bimage = zeros(sx,sy)
    bimageI = zeros(Int64,sx,sy)
    testim2 = zeros(sx,sy)
    bmaskim2 = zeros(Bool,sx,sy)
    goodpix = zeros(Bool,sx,sy)

    prelim_infill!(testim,bmaskd,bimage,bimageI,testim2,bmaskim2,goodpix;widx=19,widy=19)

    # exposure datetime based seed
    rndseed = parse(Int,date[1:6])*10^6 + parse(Int,date[8:end])
    add_sky_noise!(testim2,bmaskd,sky_im,gain;seed=rndseed)
    ## construct local covariance matrix
    # it would be nice if we handled bc well enough to not have to do the mask below
    # to do this would require cov_construct create images which are views and extend Nphalf beyond current boundaries
    # need to implement before run
    stars_interior = ((x_stars) .> Np) .& ((x_stars) .< sx.-Np) .& ((y_stars) .> Np) .& ((y_stars) .< sy.-Np);
    cxx = x_stars[stars_interior]
    cyy = y_stars[stars_interior]
    cflux = flux_stars[stars_interior]
    cov_loc, μ_loc = cov_construct(testim2, cxx, cyy; Np=Np, widx=129, widy=129)
    ## iterate over all star positions and compute errorbars/debiasing corrections
    (Nstars,) = size(cxx)
    star_stats = zeros(Nstars,7)
    for i=1:Nstars
        data_in, data_w, stars_in, kmasked2d = stamp_cutter(cxx[i],cyy[i],testim,w_im,mod_im,sky_im,bmaskd;Np=33)
        psft, kstar, kpsf2d, cntks, dnt = gen_pix_mask(kmasked2d,psfmodel,cxx[i],cyy[i],cflux[i];Np=33,thr=thr)
        star_stats[i,:] .= vec(condCovEst_wdiag(cov_loc[i,:,:],μ_loc[i,:],kstar,kpsf2d,data_in,data_w,stars_in,psft))
    end
    # prepare for export by appending to cat dictionary
    for (ind,col) in enumerate(["dcflux","dcflux_diag","dfdb","fdb","fdb_res","fdb_pred","gchi2"])
        push!(w,(col,star_stats[:,ind]))
    end

    save_fxn(Dict(w),basecat,date,filt,vers,ccd)

    return w
end

# need to write the wrapper function that loops over ccds (unfinished)
## need to think more about the memory preallocation (probably not limiting factor)
## improve commandline function access
