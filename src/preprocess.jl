module preprocess

using disCovErr
import ImageFiltering
import Distributions
import StatsBase
using Random
using LinearAlgebra

export gen_mask_staticPSF!
export prelim_infill!
export add_sky_noise!

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
        fluxt=abs(flux_stars[i])
        x_star = round(Int64, x_stars[i])
        y_star = round(Int64, y_stars[i])
        mskt = (psfstamp .> thr/fluxt)[maximum([1,2+Δx-x_star]):minimum([1+sx-x_star+Δx,psx]),maximum([1,2+Δy-y_star]):minimum([1+sy-y_star+Δy,psy])]
        bmaskd[maximum([1,x_star-Δx]):minimum([x_star+Δx,sx]),maximum([1,y_star-Δy]):minimum([y_star+Δy,sy])] .|= mskt
        # FIX ME: worth triple checking these relative indexings (remove inbounds for testing when you do that!!)
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
        disCovErr.boxsmoothMod!(bimage, in_image, wid, wid, sx, sy, 0, 0)
        disCovErr.boxsmoothMod!(bimageI, in_mask, wid, wid, sx, sy, 0, 0)

        goodpix .= (bimageI .> 10)

        testim2[bmaskim2 .& goodpix] .= (bimage./bimageI)[bmaskim2 .& goodpix]
        bmaskim2[goodpix] .= false

        # update loop params
        cnt+=1
        wid*=1.4
        wid = round(Int,wid)
        Δ = (wid-1)÷2
    end
    println("Infilling completed after $cnt rounds with final width $wid")

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

end
