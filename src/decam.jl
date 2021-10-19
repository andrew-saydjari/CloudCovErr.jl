## Handler for reading outputs of crowdsource processing on DECaPS
module decam

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

import disCovErr
import FITSIO
import ImageFiltering
import Distributions
import StatsBase
using Random
using LinearAlgebra
using PyCall
import Conda
import Pkg

"""
    __int__()

Builds the required python crowdsource dependency and exports the required load
function for obtaining the position dependent psf to the python namespace.
"""
function __init__()
    Conda.pip_interop(true)
    Conda.pip("install","crowdsourcephoto")
    Pkg.build("PyCall")
    println(Conda.PYTHONDIR)
    println(Conda.ROOTENV)
    println(joinpath(dirname(pathof(PyCall))))
    py"""
    import crowdsource.psf as psfmod
    from astropy.io import fits
    from functools import partial

    def load_psfmodel(outfn, ccd, filter, pixsz=9):
        f = fits.open(outfn)
        psfmodel = psfmod.linear_static_wing_from_record(f[ccd+"_PSF"].data[0],filter=filter)
        f.close()
        psfmodel.fitfun = partial(psfmod.fit_linear_static_wing, filter=filter, pixsz=pixsz)
        return psfmodel
    """
end

# FIX ME: Is there a world where we should be using the S7 corrected
# crowdsource imports... careful on the wt sqrt in that case though
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
function read_decam(base,date,filt,vers,ccd)
    f = FITSIO.FITS(base*date*"_ooi_"*filt*"_"*vers*".fits.fz")
    ref_im = FITSIO.read(f[ccd])
    FITSIO.close(f)
    f = FITSIO.FITS(base*date*"_oow_"*filt*"_"*vers*".fits.fz")
    w_im = FITSIO.read(f[ccd])
    FITSIO.close(f)
    f = FITSIO.FITS(base*date*"_ood_"*filt*"_"*vers*".fits.fz")
    d_im = FITSIO.read(f[ccd])
    FITSIO.close(f)
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
    f = FITSIO.FITS(base*"cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits")
    x_stars = FITSIO.read(f[ccd*"_CAT"],"x")
    y_stars = FITSIO.read(f[ccd*"_CAT"],"y")
    flux_stars = FITSIO.read(f[ccd*"_CAT"],"flux")
    decapsid = FITSIO.read(f[ccd*"_CAT"],"decapsid")
    gain = FITSIO.read_key(f[ccd*"_HDR"],"GAINCRWD")[1]
    FITSIO.close(f)

    f = FITSIO.FITS(base*"mod/c4d_"*date*"_ooi_"*filt*"_"*vers*".mod.fits")
    mod_im = FITSIO.read(f[ccd*"_MOD"])
    sky_im = FITSIO.read(f[ccd*"_SKY"])
    FITSIO.close(f)
    #switch x, y order and compensate for 0 v 1 indexing between Julia and python
    return y_stars.+1, x_stars.+1, flux_stars, decapsid, gain, mod_im, sky_im
end

"""
    gen_mask_staticPSF!(maskd, psfstamp, x_stars, y_stars, flux_stars, thr=20)

Generate a mask for an input image (which is usually an image of model residuals)
that excludes the cores of stars (which are often mismodeled). In this function,
we use a fixed PSF `psfstamp` for all sources, and adjust the masking fraction based on the
stellar flux and a threshold `thr`. A more general position dependent PSF model could be
used with a slight generalization of this function, but is likely overkill for the problem
of making a mask.

# Arguments:
- `maskd`: bool image to which mask will be added (bitwise or)
- `psfstamp`: simple 2D array of a single PSF to be used for the whole image
- `x_stars`: list of source x positions
- `y_stars`: list of source y positions
- `flux_stars`: list of source fluxes
- `thr`: threshold used for flux-dependent masking
"""
function gen_mask_staticPSF!(bmaskd, psfstamp, x_stars, y_stars, flux_stars; thr=20)
    (sx, sy) = size(bmaskd)
    (psx, psy) = size(psfstamp)
    Δx = (psx-1)÷2
    Δy = (psy-1)÷2
    Nstar = size(x_stars)[1]
    # assumes x/y_star is one indexed
    for i=1:Nstar
        fluxt=flux_stars[i]
        x_star = round(Int64, x_stars[i])
        y_star = round(Int64, y_stars[i])
        mskt = (psfstamp .> thr/abs(fluxt))[maximum([1,2+Δx-x_star]):minimum([1+sx-x_star+Δx,psx]),maximum([1,2+Δy-y_star]):minimum([1+sy-y_star+Δy,psy])]
        bmaskd[maximum([1,x_star-Δx]):minimum([x_star+Δx,sx]),maximum([1,y_star-Δy]):minimum([y_star+Δy,sy])] .|= mskt
        # FIX ME: worth triple checking these relative indexings
    end
end

"""
    prelim_infill!(testim,maskim,bimage,bimageI,testim2, maskim2, goodpix; widx = 19, widy=19)

This intial infill replaces masked pixels with a guess based on a smoothed
boxcar. For large masked regions, the smoothing scale is increased. If this
iteration takes too long/requires too strong of masking, the masked pixels
are replaced with the median of the image.

We use 3 copies of the input image and mask image. The first
is an untouched view (with reflective boundary condition padding), the second
is allocated to hold various smoothings of the image, and the third holds the
output image which contains our best infill guess. A final bool array of size
corresponding to the image is used to keep track of pixels that have safe
infill values.

# Arguments:
- `testim`: input image which requires infilling
- `bimage`: preallocated array for smoothed version of input image
- `testim2`: inplace modified ouptut array for infilled version of image input
- `maskim`: input mask indicating which pixels require infilling
- `bimageI`: preallocated array for smoothed mask counting the samples for each estimate
- `maskim2`: inplace modified mask to keep track of which pixels still need infilling
- `widx::Int`: size of boxcar smoothing window in x
- `widy::Int`: size of boxcar smoothing window in y
"""
function prelim_infill!(testim,bmaskim,bimage,bimageI,testim2,bmaskim2,goodpix;widx=19,widy=19)

    wid = maximum([widx, widy])
    Δ = (wid-1)÷2
    (sx, sy) = size(testim)

    #the masked entries in testim must be set to 0 so they drop out of the mean
    testim[bmaskim] .= 0;
    bmaskim2 .= copy(bmaskim)
    testim2 .= copy(testim)

    #loop to try masking at larger and larger smoothing to infill large holes
    cnt=0
    while any(bmaskim2) .& (cnt .< 10)
        in_image = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(Δ+2,Δ+2)));
        in_mask = ImageFiltering.padarray(.!bmaskim,ImageFiltering.Pad(:reflect,(Δ+2,Δ+2)));
        # FIX ME: could do this better to use widx != widy
        disCovErr.boxsmoothMod!(bimage, in_image, wid, wid, sx, sy)
        disCovErr.boxsmoothMod!(bimageI, in_mask, wid, wid, sx, sy)

        goodpix .= (bimageI .> 10)

        testim2[bmaskim2 .& goodpix] .= (bimage./bimageI)[bmaskim2 .& goodpix]
        bmaskim2[goodpix] .= false

        # update loop params
        cnt+=1
        wid*=1.4
        wid = round(Int,wid)
        Δ = (wid-1)÷2
    end
    println((cnt,wid))

    #catastrophic failure fallback
    if cnt == 10
        testim2[bmaskim2] .= StatsBase.median(testim)
        println("Infilling Failed Badly")
    end
    return
end

"""
    add_sky_noise!(testim2,maskim0,skyim3,gain;seed=2021)

Adds noise to the infill that matches the Poisson noise of a rough estimate for
the sky background. A random seed to set a local random generator is provided for
reproducible unit testing.

# Arguments:
- `testim2`: input image which had infilling
- `maskim0`: mask of pixels which were infilled
- `skyim3`: rough estimate of sky background counts
- `gain`: gain of detector to convert from photon count noise to detector noise
- `seed`: random seed for random generator
"""
function add_sky_noise!(testim2,maskim0,skyim3,gain;seed=2021)
    rng = MersenneTwister(seed)
    for j=1:size(testim2)[2], i=1:size(testim2)[1]
        if maskim0[i,j]
            intermed = -(rand(rng, Distributions.Poisson(convert(Float64,gain*(skyim3[i,j]-testim2[i,j]))))/gain.-skyim3[i,j])
            testim2[i,j] = intermed
        end
    end
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

"""
    stamp_cutter(cxx,cyy,residimIn,w_im,mod_im,skyim,maskim;Np=33) -> data_in, data_w, stars_in, kmasked2d

Cuts out local stamps around each star of the various input images to be used for
per star statistics calculations.

# Arguments:
- `cxx`: center coorindate x of the stamp
- `cyy`: center coorindate y of the stamp
- `residimIn`: residual image with infilling from which covariance was estimated
- `w_im`: input weight image
- `mod_im`: input model image
- `skyim`: input image of sky background
- `maskim`: input image of masked pixels
- `Np`: size of covariance matrix footprint around each star
"""
function stamp_cutter(cxx,cyy,residimIn,w_im,mod_im,skyim,maskim;Np=33)
    cx = round(Int64,cxx)
    cy = round(Int64,cyy)
    radNp = (Np-1)÷2

    cov_stamp = cx-radNp:cx+radNp,cy-radNp:cy+radNp
    @views data_in = residimIn[cov_stamp[1],cov_stamp[2]]
    @views data_w = w_im[cov_stamp[1],cov_stamp[2]]
    @views stars_in = (mod_im.-skyim)[cov_stamp[1],cov_stamp[2]]
    @views kmasked2d = maskim[cov_stamp[1],cov_stamp[2]];

    return data_in, data_w, stars_in, kmasked2d
end

function gen_pix_mask(kmasked2d,psfmodel,x_star,y_star,flux_star;Np=33,thr=thr)

    psft = psfmodel(x_star,y_star,Np)

    if flux_star < 1e4  #these are the pixels we want the cov of
        kpsf2d = (psft .> thr/1e4)
    else
        kpsf2d = (psft .> thr/flux_star)
    end

    kstar = (kmasked2d .| kpsf2d)[:]
    cntks = count(kstar)

    dnt = 0
    if cntks > 33^2-128 #this is about a 10% cut, and is the sum of bndry
        dnt = 1
        kmasked2d[1,:] .= 0
        kmasked2d[end,:] .= 0
        kmasked2d[:,1] .= 0
        kmasked2d[:,end] .= 0

        kpsf2d[1,:] .= 0
        kpsf2d[end,:] .= 0
        kpsf2d[:,1] .= 0
        kpsf2d[:,end] .= 0

        kstar = (kmasked2d .| kpsf2d)[:]
    end
    return psft, kstar, kpsf2d, cntks, dnt
end

"""
    condCovEst_wdiag(cov_loc,μ,k,kstar,kpsf2d,data_in,data_w,stars_in) -> [std_w std_wdiag var_wdb resid_mean pred_mean chi20]

Using a local covariance matrix estimate `cov_loc` and a set of known pixels `k`
and unknown pixels `kstar`, this function computes a prediction for the mean value
of the `kstar` pixels and the covariance matrix of the `kstar` pixels. In terms of
statistics use to adjust the photometry of a star, we are only interested in the
pixels masked as a result of the star (i.e. not a detector defect or cosmic ray nearby)
which is `kpsf2d`. The residual image `data_in`, the weight image `data_w`, and a model of the counts
above the background coming from the star `stars_in` for the local patch are also
inputs of the function. Correction factors for the photometric flux and flux
uncertainities are outputs as well as the chi2 value for the predicted pixels.

# Arguments:
- `cov_loc`: local covariance matrix
- `μ`: vector containing mean value for each pixel in the patch
- `k`: unmasked pixels
- `kstar`: masked pixels
- `kpsf2d`: pixels masked due to the star of interest
- `data_in`: residual image in local patch
- `data_w`: weight image in local patch
- `stars_in`: image of counts from star alone in local patch
"""
function condCovEst_wdiag(cov_loc,μ,kstar,kpsf2d,data_in,data_w,stars_in,psft)
    k = .!kstar
    kpsf1d = kpsf2d[:]
    kpsf1d_kstar = kpsf1d[kstar]
    #think about the gspice trick and invert cov_r, and then get cov_kk for free (condition on the kstar/k ratio)
    cov_r = Symmetric(cov_loc) + diagm(0 => stars_in[:])
    cov_kk = Symmetric(cov_r[k,k])
    cov_kstark = cov_r[kstar,k];
    cov_kstarkstar = Symmetric(cov_r[kstar,kstar]);
    icov_kk = inv(cholesky(cov_kk))
    predcovar = Symmetric(cov_kstarkstar - (cov_kstark*icov_kk*cov_kstark'))
    ipcov = inv(cholesky(predcovar))

    @views uncond_input = data_in[:]
    @views cond_input = data_in[:].- μ

    kstarpredn = cov_kstark*icov_kk*cond_input[k]
    kstarpred = kstarpredn .+ μ[kstar]
    #This is not that impressive. chi2 of larger patch is more useful diagnostic
    chi20 = (kstarpredn'*ipcov*kstarpredn)

    @views p = psft[kpsf2d][:]
    @views pw = p.*data_w[kpsf2d][:]
    @views p2w = p.*pw

    @views std_w = sqrt(abs(pw'*predcovar[kpsf1d_kstar,kpsf1d_kstar]*pw))/sum(p2w)
    @views std_wdiag = sqrt(abs(sum((pw.^(2)).*diag(predcovar[kpsf1d_kstar,kpsf1d_kstar]))))/sum(p2w)
    @views var_wdb = (p'*ipcov[kpsf1d_kstar,kpsf1d_kstar]*p)

    @views resid_mean = (p'*ipcov[kpsf1d_kstar,kpsf1d_kstar]*uncond_input[kpsf1d])./var_wdb
    @views pred_mean = (p'*ipcov[kpsf1d_kstar,kpsf1d_kstar]*kstarpred[kpsf1d_kstar])./var_wdb

    return [std_w std_wdiag var_wdb resid_mean pred_mean chi20]
end
end
