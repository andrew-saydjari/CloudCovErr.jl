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
    # if !haskey(Conda._installed_packages_dict(),"crowdsourcephoto")
    #     Conda.add("nomkl",channel="conda-forge")
    #     Conda.add("crowdsourcephoto",channel="conda-forge")
    # end
    py"""
    import sys
    sys.path.append('/n/home12/saydjari/finksage/ExtSoftware/crowdsource_new/crowdsource')
    import psf as psfmod
    from decam_proc import read_data
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
    dfn = base*date*"_ood_"*filt*"_"*vers*".fits.fz"
    if last(ccd,1) == "I"
        ifn = inject_rename(ifn)
        dfn = inject_rename(dfn)
    end
    if corrects7 .& ((ccd == "S7") .| (ccd == "S7I"))
        # a little wasteful to throw away load times in python for w_im and d_im
        # but we need the crowdsource formatting for s7 correction and it is easiest
        # to keep the disCovErr formatting consistent
        wfn = base*date*"_oow_"*filt*"_"*vers*".fits.fz"
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
    f = FITS(dfn)
    d_im = read(f[ccd])
    close(f)
    return ref_im, d_im
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
    for i=1:length(wcol)
        if (wcol[i] == "passno") | (wcol[i] == "dnt")
            w[i] = convert.(Int8,w[i])
        elseif (wcol[i] == "kcond0") | (wcol[i] == "kcond") | (wcol[i] == "kpred")
            w[i] = convert.(Int32,w[i])
        elseif (wcol[i] == "rchi2") | (wcol[i] == "spread_model") | (wcol[i] == "dspread_model")
            w[i] = convert.(Float32,w[i])
        elseif wcol[i] == "fdb_tot"
            for j=1:length(wcol)
                if wcol[j] == "flux"
                    w[i].+=w[j]
                    wcol[i] = "cflux"
                end
            end
        end
    end
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

function proc_ccd(base,date,filt,vers,basecat,ccd;thr=20,Np=33,corrects7=true,widx=129,widy=widx,tilex=1,tiley=tilex,ftype::Int=32)
    println("Started $ccd")
    flush(stdout)

    if ftype == 32
        T = Float32
    else
        T = Float64
    end

    # loads from disk
    ref_im, d_im = read_decam(base,date,filt,vers,ccd,corrects7=corrects7)
    bmaskd = (d_im .!= 0)
    d_im = nothing
    (sx0, sy0) = size(ref_im)
    x_stars, y_stars, flux_stars, decapsid, gain, mod_im, sky_im, wcol, w = cloudCovErr.read_crowdsource(basecat,date,filt,vers,ccd)

    psfmodel = cloudCovErr.load_psfmodel_cs(basecat,date,filt,vers,ccd)
    psfstatic511 = psfmodel(sx0÷2,sy0÷2,511)
    psfstatic33 = psfmodel(sx0÷2,sy0÷2,Np)

    # mask bad camera pixels/cosmic rays, then mask out star centers
    cloudCovErr.gen_mask_staticPSF2!(bmaskd, psfstatic511, psfstatic33, x_stars, y_stars, flux_stars; thr=thr)

    testim = copy(mod_im .- ref_im)
    bimage = zeros(T,sx0,sy0)
    bimageI = zeros(Int64,sx0,sy0)
    testim2 = zeros(T,sx0,sy0)
    bmaskim2 = zeros(Bool,sx0,sy0)
    goodpix = zeros(Bool,sx0,sy0)

    prelim_infill!(testim,bmaskd,bimage,bimageI,testim2,bmaskim2,goodpix,ccd;widx=19,widy=19,ftype=ftype)
    testim = copy(mod_im .- ref_im) ##FIXME should not be overwritten
    ref_im = nothing
    bimage = nothing
    bimageI = nothing
    bmaskim2 = nothing
    goodpix = nothing
    return testim2
    ## calculate the star farthest outside the edge of the image in x and y
    cx = round.(Int,x_stars)
    cy = round.(Int,y_stars)
    px0 = outest_bounds(cx,sx0)
    py0 = outest_bounds(cy,sy0)

    ## these have to be allocating to get the noise model right
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    padx = Np+Δx+px0
    pady = Np+Δy+py0
    in_image = ImageFiltering.padarray(testim2,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    testim2 = nothing
    in_image_raw = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    testim = nothing
    in_sky_im = ImageFiltering.padarray(sky_im,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    in_stars_im = ImageFiltering.padarray(abs.(mod_im.-sky_im)./gain,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    sky_im = nothing
    mod_im = nothing
    in_bmaskd = ImageFiltering.padarray(bmaskd,ImageFiltering.Fill(true,(padx+2,pady+2)));
    bmaskd = nothing

    # exposure datetime based seed
    rndseed = parse(Int,date[1:6])*10^6 + parse(Int,date[8:end])
    cloudCovErr.add_sky_noise!(in_image,in_bmaskd,in_sky_im,gain;seed=rndseed)
    in_sky_im = nothing

    ## iterate over all star positions and compute errorbars/debiasing corrections
    (Nstars,) = size(x_stars)
    star_stats = zeros(T,10,Nstars)

    # preallocate the cov and μ per star variables
    cov = zeros(T,Np*Np,Np*Np)
    μ = zeros(T,Np*Np)

    # compute a radial mask for reduced num cond pixels
    circmask = kstar_circle_mask(Np,rlim=256)

    # some important global sizes for the loop
    cntStar0 = 0
    stepx = (sx0+2) ÷ tilex
    stepy = (sy0+2) ÷ tiley

    # precallocate the image subblocks
    in_subimage = zeros(T,stepx+2*padx,stepy+2*pady)
    ism = zeros(T,stepx+2*padx,stepy+2*pady)
    bimage = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy)
    bism = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,2*Np-1, Np);
    for jx=1:tilex, jy=1:tiley
        xrng, yrng, star_ind = im_subrng(jx,jy,cx,cy,sx0+2,sy0+2,px0,py0,stepx,stepy,padx,pady,tilex,tiley)
        cntStar = length(star_ind)
        if cntStar > 0
            in_subimage .= in_image[xrng,yrng]
            cov_avg!(bimage, ism, bism, in_subimage, widx=widx, widy=widy)
            offx = padx-Δx-(jx-1)*stepx
            offy = pady-Δy-(jy-1)*stepy
            for i in star_ind
                build_cov!(cov,μ,cx[i]+offx,cy[i]+offy,bimage,bism,Np,widx,widy)
                data_in, stars_in, kmasked2d = stamp_cutter(cx[i],cy[i],in_image_raw,in_stars_im,in_bmaskd;Np=Np)
                psft, kstar, kpsf2d, kcond0, kcond, kpred, dnt = gen_pix_mask(kmasked2d,psfmodel,circmask,x_stars[i],y_stars[i],flux_stars[i];Np=Np,thr=thr)
                try
                    star_stats[:,i] .= [condCovEst_wdiag(cov,μ,kstar,kpsf2d,data_in,stars_in,psft,diag_on=true)[1]..., kcond0, kcond, kpred, dnt]
                catch
                    star_stats[:,i] .= [NaN, NaN, NaN, NaN, NaN, NaN, kcond0, kcond, kpred, dnt]
                end
            end
        end
        cntStar0 += cntStar
        println("Finished $cntStar stars in tile ($jx, $jy) of $ccd")
        flush(stdout)
    end
    # prepare for export by appending to cat vectors
    # if doing the ops float64, might want to do a final convert to float32 before saving
    for (ind,col) in enumerate(["dcflux","dcflux_diag","fdb_tot","fdb_res","fdb_pred","cchi2","kcond0","kcond","kpred","dnt"])
        push!(wcol,col)
        push!(w,star_stats[ind,:])
    end
    cloudCovErr.save_fxn(wcol,w,basecat,date,filt,vers,ccd)
    println("Saved $ccd processing $cntStar0 of $Nstars stars")
    flush(stdout)
    pdefer = count(isnan.(star_stats[1,:]))
    if pdefer > 0
        println("There were $pdefer posDef errors")
        flush(stdout)
    end

    return
end

function proc_all(base,date,filt,vers,basecat;ccdlist=String[],resume=false,corrects7=true,thr=20,Np=33,widx=129,widy=widx,tilex=1,tiley=tilex,ftype::Int=32)
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
        proc_ccd(base,date,filt,vers,basecat,ccd;thr=thr,Np=Np,corrects7=corrects7,widx=widx,widy=widx,tilex=tilex,tiley=tilex,ftype=ftype)
        GC.gc()
    end
end

end
