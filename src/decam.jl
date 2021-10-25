## Handler for reading outputs of crowdsource processing on DECaPS
module decam

export inject_rename
export read_decam
export read_crowdsource
export load_psfmodel_cs
export save_fxn
export get_catnames

export proc_ccd
export proc_all

using cloudCovErr
using PyCall
using FITSIO
import ImageFiltering
import Distributions
import StatsBase
using Random
using LinearAlgebra
import Conda

"""
    __int__()

Builds the required python crowdsource dependency and exports the required load
function for obtaining the position dependent psf to the python namespace.
"""
function __init__()
    if !haskey(Conda._installed_packages_dict(),"crowdsourcephoto")
        Conda.add("nomkl",channel="conda-forge")
        Conda.add("crowdsourcephoto",channel="conda-forge")
    end
    py"""
    import crowdsource.psf as psfmod
    from crowdsource.decam_proc import read_data
    from astropy.io import fits
    from functools import partial
    import os

    # default decam_dir at Harvard
    if 'DECAM_DIR' not in os.environ:
        os.environ['DECAM_DIR'] = '/n/home13/schlafly/decam'

    def load_psfmodel(outfn, ccd, filter, pixsz=9):
        f = fits.open(outfn)
        psfmodel = psfmod.linear_static_wing_from_record(f[ccd+"_PSF"].data[0],filter=filter)
        f.close()
        psfmodel.fitfun = partial(psfmod.fit_linear_static_wing, filter=filter, pixsz=pixsz)
        return psfmodel
    """
end

## File specific read/write functions

function inject_rename(fname)
    splitname = split(fname,"/")
    splitname[4] = "decapsi"
    return chop(join(splitname,"/"),tail=7)*"I.fits.fz"
end

"""
    read_decam(base,date,filt,vers,ccd) -> ref_im, w_im, d_im

Read in raw image files associated with exposures obtain on the DarkEnergyCamera.
Returns the image, a weighting image, and a quality flag mask image. See [NOAO
handbook](http://ast.noao.edu/sites/default/files/NOAO_DHB_v2.2.pdf) for more details
on what is contained in each file and how they are obtained.

# Arguments:
- `base`: parent directory and file name prefix for exposure files
- `date`: date_time of the exposure
- `filt`: optical filter used to take the exposure
- `vers`: NOAO community processing version number
- `ccd`: which ccd we are pulling the image for

# Example
```julia
ref_im, w_im, d_im = read_decam("/n/fink2/decaps/c4d_","170420_040428","g","v1","N14")
```

"""
function read_decam(base,date,filt,vers,ccd; corrects7=true)
    ifn = base*date*"_ooi_"*filt*"_"*vers*".fits.fz"
    wfn = base*date*"_oow_"*filt*"_"*vers*".fits.fz"
    dfn = base*date*"_ood_"*filt*"_"*vers*".fits.fz"
    if last(ccd,1) == "I"
        ifn = inject_rename(ifn)
        wfn = inject_rename(wfn)
        dfn = inject_rename(dfn)
    end
    if corrects7 .& ((ccd == "S7") .| (ccd == "S7I"))
        # a little wasteful to throw away load times in python for w_im and d_im
        # but we need the crowdsource formatting for s7 correction and it is easiest
        # to keep the disCovErr formatting consistent
        py_ref_im, py_w_im, py_d_im, nebprob = py"read_data"(ifn,wfn,dfn,ccd,maskdiffuse=false)
        py_w_im = nothing
        py_d_im = nothing
        nebprob = nothing
        ref_im = py_ref_im'
    else
        f = FITS(ifn)
        ref_im = read(f[ccd])
        close(f)
    end
    f = FITS(wfn)
    w_im = read(f[ccd])
    close(f)
    f = FITS(dfn)
    d_im = read(f[ccd])
    close(f)
    return ref_im, w_im, d_im
end

"""
    read_crowdsource(base,date,filt,vers,ccd) -> x_stars, y_stars, flux_stars, decapsid, gain, mod_im, sky_im

Read in outputs of crowdsource, a photometric pipeline. To pair with an arbitrary
photometric pipeline, an analogous read in function should be created. The relevant
outputs are the model image (including the sources) so that we can produce the
residual image, the sky/background model (no sources), and the coordinates of the stars.
The survey id number is also readout of the pipeline solution file to help
cross-validate matching of the disCovErr outputs and the original sources. The empirical
gain is read out of the header (for other photometric pipelines which don't perform this estiamte,
the gain from DECam is likely sufficient).

# Arguments:
- `base`: parent directory and file name prefix for crowdsource results files
- `date`: date_time of the exposure
- `filt`: optical filter used to take the exposure
- `vers`: NOAO community processing version number
- `ccd`: which ccd we are pulling the image for
"""
function read_crowdsource(base,date,filt,vers,ccd)
    f = FITS(base*"cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits")
    x_stars = read(f[ccd*"_CAT"],"x")
    y_stars = read(f[ccd*"_CAT"],"y")
    flux_stars = read(f[ccd*"_CAT"],"flux")
    decapsid = read(f[ccd*"_CAT"],"decapsid")
    gain = read_key(f[ccd*"_HDR"],"GAINCRWD")[1]
    wcol = String[]
    w = []
    for col in FITSIO.colnames(f[ccd*"_CAT"])
        push!(wcol,col)
        push!(w,read(f[ccd*"_CAT"],col))
    end
    close(f)

    f = FITS(base*"mod/c4d_"*date*"_ooi_"*filt*"_"*vers*".mod.fits")
    mod_im = read(f[ccd*"_MOD"])
    sky_im = read(f[ccd*"_SKY"])
    close(f)
    #switch x, y order and compensate for 0 v 1 indexing between Julia and python
    return y_stars.+1, x_stars.+1, flux_stars, decapsid, gain, mod_im, sky_im, wcol, w
end

"""
    load_psfmodel_cs(base,date,filt,vers,ccd) -> psfmodel

Julia wrapper function for the PyCall that reads the position dependent psfmodel
produced by crowdsource from the catalogue file for a given exposure and ccd. The
returned psfmodel takes an x- and y-position for the source location and the size
of the desired psfstamp (the stamps are square and required to be odd).

# Arguments:
- `base`: parent directory and file name prefix for catalogue files
- `date`: date_time of the exposure
- `filt`: optical filter used to take the exposure
- `vers`: NOAO community processing version number
- `ccd`: which ccd we are pulling the image for
"""
function load_psfmodel_cs(base,date,filt,vers,ccd)
    psfmodel_py = py"load_psfmodel"(base*"cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits",ccd,filt)
    function psfmodel_jl(x,y,sz)
        # accounts for x, y ordering and 0 v 1 indexing between python and Julia
        return psfmodel_py(y.-1,x.-1,sz)'
    end
    return psfmodel_jl
end

function save_fxn(wcol,w,base,date,filt,vers,ccd)
    f = FITS(base*"cer/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.cer.fits","r+")
    write(f,wcol,w,name=ccd*"_CAT")
    close(f)
end

function get_catnames(f)
    extnames = String[]
    i=1
    for h in f
        if i==1
        else
            extname = read_key(h,"EXTNAME")[1]
            if last(extname,3) == "CAT"
                push!(extnames,chop(extname,tail=4))
            end
        end
        i+=1
    end
    return extnames
end

## Main run functions for processing exposures
# combines all of the functions in the repo for
# an acutal implementation

# should we be prellocating outside this subfunction?
function proc_ccd(base,date,filt,vers,basecat,ccd;thr=20,Np=33,corrects7=true)
    println("Started $ccd")
    # loads from disk
    ref_im, w_im, d_im = read_decam(base,date,filt,vers,ccd,corrects7=corrects7)
    (sx, sy) = size(ref_im)
    x_stars, y_stars, flux_stars, decapsid, gain, mod_im, sky_im, wcol, w = read_crowdsource(basecat,date,filt,vers,ccd)

    psfmodel = load_psfmodel_cs(basecat,date,filt,vers,ccd)
    psfstatic = psfmodel(sx÷2,sy÷2,511)

    # mask bad camera pixels/cosmic rays, then mask out star centers
    bmaskd = (d_im .!= 0)
    gen_mask_staticPSF!(bmaskd, psfstatic, x_stars, y_stars, flux_stars; thr=thr)

    testim = copy(mod_im .- ref_im)
    bimage = zeros(Float64,sx,sy)
    bimageI = zeros(Int64,sx,sy)
    testim2 = zeros(Float64,sx,sy)
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
    #stars_interior = ((x_stars) .> Np) .& ((x_stars) .< sx.-Np) .& ((y_stars) .> Np) .& ((y_stars) .< sy.-Np);
    cxx = x_stars #[stars_interior]
    cyy = y_stars #[stars_interior]
    cflux = flux_stars #[stars_interior]
    cov_loc, μ_loc = cov_construct(testim2, cxx, cyy; Np=Np, widx=129, widy=129)
    ## iterate over all star positions and compute errorbars/debiasing corrections
    (Nstars,) = size(cxx)
    star_stats = zeros(Nstars,7)

    in_image = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(Np+2,Np+2)));
    in_w_im = ImageFiltering.padarray(w_im,ImageFiltering.Pad(:reflect,(Np+2,Np+2)));
    in_mod_im = ImageFiltering.padarray(mod_im,ImageFiltering.Pad(:reflect,(Np+2,Np+2)));
    in_sky_im = ImageFiltering.padarray(sky_im,ImageFiltering.Pad(:reflect,(Np+2,Np+2)));
    in_bmaskd = ImageFiltering.padarray(bmaskd,ImageFiltering.Pad(:reflect,(Np+2,Np+2)));

    for i=1:Nstars
        data_in, data_w, stars_in, kmasked2d = stamp_cutter(cxx[i],cyy[i],in_image,in_w_im,in_mod_im,in_sky_im,in_bmaskd;Np=33)
        psft, kstar, kpsf2d, cntks, dnt = gen_pix_mask(kmasked2d,psfmodel,cxx[i],cyy[i],cflux[i];Np=33,thr=thr)
        star_stats[i,:] .= vec(condCovEst_wdiag(cov_loc[i,:,:],μ_loc[i,:],kstar,kpsf2d,data_in,data_w,stars_in,psft))
    end
    # prepare for export by appending to cat vectors
    for (ind,col) in enumerate(["dcflux","dcflux_diag","dfdb","fdb","fdb_res","fdb_pred","gchi2"])
        push!(wcol,col)
        push!(w,star_stats[:,ind])
    end
    save_fxn(wcol,w,basecat,date,filt,vers,ccd)
    println("Saved $ccd")
    return
end

function proc_all(base,date,filt,vers,basecat;ccdlist=String[],resume=false,corrects7=true,thr=20,Np=33)
    infn = basecat*"cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits"
    println("Starting to process "*infn)

    f = FITS(infn)
    prihdr = read_header(f[1])
    extnames=get_catnames(f)
    close(f)

    outfn = basecat*"cer/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.cer.fits"

    # write output file or read it in and check which ccds are completed
    if ((!resume) | (!isfile(outfn)))
        f = FITS(outfn,"w")
        write(f,[0], header=prihdr)
        close(f)
        extnamesdone = String[]
    else
        f = FITS(outfn,"r+")
        extnamesdone=get_catnames(f)
        close(f)
    end

    ## prepare list of ccds to run
    # remove the ones already done
    if length(extnamesdone) !=0
        alreadydone = intersect(extnames,extnamesdone)
        extnames = setdiff(extnames,extnamesdone)
        println("Skipping ccds already completed: ", alreadydone)
    end
    # limit to the ccds in the list not yet run
    if length(ccdlist) != 0
        extnames = intersect(extnames,ccdlist)
        println("Only running unfinished ccds in ccdlist: ", extnames)
    end

    # main loop over ccds
    for ccd in extnames
        proc_ccd(base,date,filt,vers,basecat,ccd;thr=thr,Np=Np,corrects7=corrects7)
    end
end

end
