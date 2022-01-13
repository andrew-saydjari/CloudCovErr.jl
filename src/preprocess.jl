using cloudCovErr
import ImageFiltering
import Distributions
import StatsBase
using Random
using LinearAlgebra

export gen_mask_staticPSF!
export gen_mask_staticPSF2!
export prelim_infill!
export add_sky_noise!
export findmaxpsf
export kstar_circle_mask
export im_subrng

"""
    kstar_circle_mask(Np;rlim=256) -> circmask

Generates a Bool mask for pixels beyond a given (squared) radius of the center of an image.

# Arguments:
- `Np`: size of image stamp

# Keywords:
- `rlim`: squared radius (in pixels^2) beyond which pixels should be masked (default 256)

# Output:
- `circmask`: static Bool mask used for assigning pixels beyond some radius of the stellar center as "ignored"
"""
function kstar_circle_mask(Np;rlim=256)
    halfNp = (Np-1) ÷ 2
    x = (-halfNp:halfNp)' .* ones(Int,Np)
    y = ones(Int,Np)' .* (-halfNp:halfNp)
    R = x.^2 .+ y.^2
    return R.>rlim
end

"""
    findmaxpsf(psfstamp1;thr=20) -> flim

Computes the flux a star must have so that the PSF-based masking using `thr`
would require a larger stamp area. Used for computational savings.

# Arguments:
- `psfstamp1`: a small image of a representative PSF

# Keywords:
- `thr`: threshold used to determine which pixels are bright enough to be "hidden"

# Output:
- `flim`: maximum flux that can be masked by `thr` without exceeding the PSF stamp footprint
"""
function findmaxpsf(psfstamp1;thr=20)
    thr/maximum([psfstamp1[:,1]...,psfstamp1[1,:]...,psfstamp1[:,end]...,psfstamp1[end,:]...])
end

"""
    im_subrng(jx,jy,cx,cy,sx,sy,px0,py0,stepx,stepy,padx,pady,tilex,tiley) -> xrng, yrng, star_ind

Computes the flux a star must have so that the PSF-based masking using `thr`
would require a larger stamp area. Used for computational savings.

# Arguments:
- `jx`: tile index along x
- `jy`: tile index along y
- `cx`: list of stellar x-coordinates
- `cy`: list of stellar y-coordinates
- `sx`: size of original image in x
- `sy`: size of original image in y
- `px0`: maximal padding in x to account for stars outside image
- `py0`: maximal padding in y to account for stars outside image
- `stepx`: tiling step size in x
- `stepy`: tiling step size in y
- `padx`: tile padding in x required to account for local stamp size, sample size, and pixels outside the image
- `pady`: tile padding in x required to account for local stamp size, sample size, and pixels outside the image
- `tilex`: total number of tile divisions along x
- `tiley`: total number of tile divisions along y

# Output:
- `xrng`: slicing range of the tile in x
- `yrng`: slicing range of the tile in y
- `star_ind`: Bool mask of all stars falling within the tile (subimage)
"""
function im_subrng(jx,jy,cx,cy,sx,sy,px0,py0,stepx,stepy,padx,pady,tilex,tiley)
    lowbx = (1 + (jx-1)*stepx)
    uppbx = (1 + jx*stepx-1)
    lowby = (1 + (jy-1)*stepy)
    uppby = (1 + jy*stepy-1)
    # pad right and top by 2 to fix ragged
    if uppbx > sx
        uppbx = sx
    end
    if uppby > sy
        uppby = sy
    end
    xrng = (lowbx-padx):(uppbx+padx)
    yrng = (lowby-pady):(uppby+pady)

    if jx == 1
        lowbx-=px0
    end
    if jx == tilex
        uppbx+=px0
    end
    if jy == 1
        lowby-=py0
    end
    if jy == tiley
        uppby+=py0
    end

    star_ind = findall((lowbx-0.5 .< cx .<= uppbx+0.5) .& (lowby-0.5 .< cy .<= uppby+0.5))

    return xrng, yrng, star_ind
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

# Keywords:
- `thr`: threshold used for flux-dependent masking (default 20)
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
        @views mskt = (psfstamp .> thr/fluxt)[maximum([1,2+Δx-x_star]):minimum([1+sx-x_star+Δx,psx]),maximum([1,2+Δy-y_star]):minimum([1+sy-y_star+Δy,psy])]
        @views bmaskd[maximum([1,x_star-Δx]):minimum([x_star+Δx,sx]),maximum([1,y_star-Δy]):minimum([y_star+Δy,sy])] .|= mskt
    end
end

"""
    gen_mask_staticPSF2!(maskd, psfstamp, psfstamp1, x_stars, y_stars, flux_stars, thr=20)

Generate a mask for an input image (which is usually an image of model residuals)
that excludes the cores of stars (which are often mismodeled). In this function,
we use a small fixed PSF `psfstamp1` for all faint sources, and adjust the masking
fraction based on the stellar flux and a threshold `thr`. Only for source bright
enough to need a larger PSF stamp do we use `psfstamp`, which saves some computational
cost.

# Arguments:
- `maskd`: bool image to which mask will be added (bitwise or)
- `psfstamp`: simple 2D array of a single PSF to be used for bright stars in the whole image
- `psfstamp1`: simple 2D array of a single PSF to be used for faint stars in the whole image
- `x_stars`: list of source x positions
- `y_stars`: list of source y positions
- `flux_stars`: list of source fluxes

# Keywords:
- `thr`: threshold used for flux-dependent masking (default 20)
"""
function gen_mask_staticPSF2!(bmaskd, psfstamp, psfstamp1, x_stars, y_stars, flux_stars; thr=20)
    (sx, sy) = size(bmaskd)
    (psx, psy) = size(psfstamp)
    Δx = (psx-1)÷2
    Δy = (psy-1)÷2
    (psx1, psy1) = size(psfstamp1)
    Δx1 = (psx1-1)÷2
    Δy1 = (psy1-1)÷2
    Nstar = size(x_stars)[1]
    flim = findmaxpsf(psfstamp1;thr=thr)
    # assumes x/y_star is one indexed
    for i=1:Nstar
        fluxt = abs(flux_stars[i])
        x_star = round(Int64, x_stars[i])
        y_star = round(Int64, y_stars[i])
        if fluxt > flim
            pxrange = maximum([1,2+Δx-x_star]):minimum([1+sx-x_star+Δx,psx])
            pyrange = maximum([1,2+Δy-y_star]):minimum([1+sy-y_star+Δy,psy])
            bxrange = maximum([1,x_star-Δx]):minimum([x_star+Δx,sx])
            byrange = maximum([1,y_star-Δy]):minimum([y_star+Δy,sy])
            @views mskt = (psfstamp[pxrange,pyrange] .> thr/fluxt)
        else
            pxrange = maximum([1,2+Δx1-x_star]):minimum([1+sx-x_star+Δx1,psx1])
            pyrange = maximum([1,2+Δy1-y_star]):minimum([1+sy-y_star+Δy1,psy1])
            bxrange = maximum([1,x_star-Δx1]):minimum([x_star+Δx1,sx])
            byrange = maximum([1,y_star-Δy1]):minimum([y_star+Δy1,sy])
            @views mskt = (psfstamp1[pxrange,pyrange] .> thr/fluxt)
        end
        #println((size(mskt),size(bmaskd[bxrange,byrange])))
        @views bmaskd[bxrange,byrange] .|= mskt
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
function prelim_infill!(testim,bmaskim,bimage,bimageI,testim2,bmaskim2,goodpix,ccd;widx=19,widy=19,ftype::Int=32,widmult=1.4)
    if ftype == 32
        T = Float32
    else
        T = Float64
    end
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    (sx, sy) = size(testim)

    widxMax = round(Int,(widmult^10)*(widx-1)/2)*2+1
    widyMax = round(Int,(widmult^10)*(widy-1)/2)*2+1
    ΔxMax = (widxMax-1)÷2
    ΔyMax = (widyMax-1)÷2

    #the masked entries in testim must be set to 0 so they drop out of the mean
    testim[bmaskim] .= 0;
    bmaskim2 .= copy(bmaskim)
    testim2 .= copy(testim)

    #hopefully replace with the reflected indexedx arrays
    in_image = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(ΔxMax,ΔxMax)));
    in_mask = ImageFiltering.padarray(.!bmaskim,ImageFiltering.Pad(:reflect,(ΔyMax,ΔyMax)));

    #loop to try masking at larger and larger smoothing to infill large holes
    cnt=0
    while any(bmaskim2) .& (cnt .< 10)
        # this double deep view should be only 1 deep ideally... need the internal unwrap
        @views in_image1 = in_image[(1-Δx):(sx+Δx),(1-Δy):(sy+Δy)]
        @views in_mask1 = in_mask[(1-Δx):(sx+Δx),(1-Δy):(sy+Δy)]
        (sx1, sy1) = size(in_image1)
        tot = zeros(T,sx1)
        totI = zeros(Int,sx1)
        boxsmooth!(bimage,in_image1,tot,widx,widy)
        boxsmooth!(bimageI,in_mask1,totI,widx,widy)

        goodpix .= (bimageI .> 10)

        testim2[bmaskim2 .& goodpix] .= (bimage./bimageI)[bmaskim2 .& goodpix]
        bmaskim2[goodpix] .= false

        # update loop params
        cnt+=1
        widx*=1.4
        widy*=1.4
        widx = round(Int,(widx-1)/2)*2+1
        widy = round(Int,(widy-1)/2)*2+1
        Δx = (widx-1)÷2
        Δy = (widy-1)÷2
    end
    println("Infilling $ccd completed after $cnt rounds with final width (widx,widy) = ($widx,$widy)")
    flush(stdout)

    #catastrophic failure fallback
    if cnt == 10
        testim2[bmaskim2] .= StatsBase.median(testim)
        println("Infilling Failed Badly")
        flush(stdout)
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
function add_sky_noise!(testim2,maskim,skyim,gain;seed=2021)
    rng = MersenneTwister(seed)
    for i in eachindex(testim2)
        if maskim[i]
            intermed = -(rand(rng, Distributions.Poisson(convert(Float64,gain*(abs.(skyim[i]-testim2[i])))))/gain.-skyim[i])
            testim2[i] = intermed
        end
    end
end

function add_noise!(testim2,gain;seed=2021)
    rng = MersenneTwister(seed)
    for i in eachindex(testim2)
        @inbounds testim2[i] = rand(rng, Distributions.Poisson(convert(Float64,gain*testim2[i])))/gain
    end
end
