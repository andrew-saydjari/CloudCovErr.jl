## utility functions
import OffsetArrays
import ImageFiltering
import ShiftedArrays
import OffsetArrays

export cov_construct
export boxsmoothMod!

"""
    cov_construct(img, cx, cy; Np::Int=33, widx::Int=129, widy::Int=129) -> cov, μ

Construct a local covariance matrix estimate from `img` centered around pixel `cx` and
`cy`. The covariance matrix will be for a square subimage of size `Np` by `Np` pixels,
yielding `cov` of size `Np`^2. The covariance matrix is estimated by taking samples over
pixels within a box of size `widx` and `widy` centered on the pixel of interest.

If `cx` and `cy` are arrays, the first index of the returned `cov` will be the index
of the center point and `μ` will be the same length as `cx`.

In the interest of speed for the case of wanting to construct the covariance matrix
at many posiitons, the moving mean of the full input image (times shifts of the input
image) is precomputed for every entry of the covariance matrix. This may not be efficient
for few `cx`.

Note that the values of `widx` and `widy` determine how local an estimate the covariance matrix
returned is, but for too small of values of `widx` and `widy`, the covariance matrix estimate
may be singular (one should always keep `Np`^2 < `widx`*`widy`).

# Arguments:
- `img`: Input image for which we desire the covariance matrix estimate
- `cx`: x position (or positions) at which we want a covariance matrix estimate
- `cy`: y position (or positions) at which we want a covariance matrix estimate
- `Np`: size of the image stamp for which the covariance matrix will be estimated
- `widx`: size of sampled region in x used to estimate entries in the covariance matrix
- `widy`: size of sampled region in y used to estimate the entries in the covariance matrix
"""
function cov_construct(img, cxx, cyy; Np::Int=33, widx::Int=129, widy::Int=129)
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    halfNp = (Np-1) ÷ 2
    (Nstar,) = size(cxx)
    (sx, sy) = size(img)

    cx = round.(Int,cxx)
    cy = round.(Int,cyy)

    px0 = maximum(abs.(cx))
    py0 = maximum(abs.(cy))

    # preallocate output covariance array
    Δr = zeros(Int,Nstar)
    Δc = zeros(Int,Nstar)
    # FIX ME: do we really benefit from this being full Float64
    cov = zeros(Nstar,Np*Np,Np*Np)
    μ = zeros(Nstar,Np*Np)

    in_image = ImageFiltering.padarray(img,ImageFiltering.Pad(:reflect,(Np+Δx+2,Np+Δy+2)));
    bism = ImageFiltering.padarray(copy(img),ImageFiltering.Pad(:reflect,(halfNp+px0,halfNp+py0)));
    bimage = ImageFiltering.padarray(copy(img),ImageFiltering.Pad(:reflect,(halfNp+px0,halfNp+py0)));
    #this padding could be smaller I think for the b_ series

    Δr, Δc = cx.-(halfNp-1), cy.-(halfNp-1)

    boxsmoothMod!(bimage,in_image,widx,widy,sx,sy,halfNp+px0,halfNp+py0)
    # loop over shifts
    for dc=0:Np-1       # column shift loop
        pcr = 1:Np-dc
        for dr=1-Np:Np-1   # row loop, incl negatives
            if (dr < 0) & (dc == 0)
                continue
            end
            if dr >= 0
                prr = 1:Np-dr
            end
            if (dr < 0) & (dc > 0)
                prr = 1-dr:Np
            end
            # ism = image, shifted and multipled
            ism = in_image .* OffsetArrays.OffsetArray(ShiftedArrays.circshift(in_image.parent,(-dr, -dc)), OffsetArrays.Origin(in_image.offsets.+1))
            boxsmoothMod!(bism,ism,widx,widy,sx,sy,halfNp+px0,halfNp+py0) # bism = boxcar(ism)

            for pc=pcr, pr=prr
                i = ((pc   -1)*Np)+pr
                j = ((pc+dc-1)*Np)+pr+dr
                for st=1:Nstar
                    drr = Δr[st]
                    dcc = Δc[st]
                    μ1μ2 = bimage[pr+drr,pc+dcc]*bimage[pr+dr+drr,pc+dc+dcc]
                    cov[st,i,j] = bism[pr+drr,pc+dcc]/(widx*widy) - μ1μ2/((widx*widy)^2)
                    if i == j
                        μ[st,i] = μ1μ2/((widx*widy)^2)
                    end
                end
            end

        end
    end
    return (widx*widy)/((widx*widy)-1).*cov, μ # before using, must call Symmetric(cov[i,:,:])
end

"""
    boxsmoothMod!(out, arr, widx::Int, widy::Int, sx::Int, sy::Int)

Boxcar smooths an input image (or paddedview) `arr` with window size `widx` by
`widy`. We pass the original image size `sx` and `sy` to help handle image views.

# Arguments:
- `out`: preallocated output array for the boxcar smoothed image
- `arr`: input array for which boxcar smoothing is computed (generally paddedview)
- `widx::Int`: size of boxcar smoothing window in x
- `widy::Int`: size of boxcar smoothing window in y
- `sx::Int`: x size of the original (unpadded) image
- `sy::Int`: y size of the original (unpadded) image
"""
function boxsmoothMod!(out, arr, widx::Int, widy::Int, sx::Int, sy::Int, px::Int, py::Int)
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2

    tot = zeros(sy)

    @inbounds for j=(1-Δy-py):(sy-Δy+py)
        if (j==1-Δy-py)
            tot = (sum(arr[:,1-Δy-py:1+Δy-py], dims=2))[:,1]
        else
            tot .+= (arr[:,j+widy-1]-arr[:,j-1])
        end
        tt=0
        @inbounds for i=(1-Δx-px):(sx-Δx+px)
            if (i==1-Δx-px)
                tt = sum(tot[1-Δx-px:1+Δx-px])
            else
                tt -= (tot[i-1]-tot[i+widx-1])
            end
            out[i+Δx,j+Δy] = tt
        end
    end
end
