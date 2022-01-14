## utility functions
import OffsetArrays
import ImageFiltering
import ShiftedArrays

export cov_avg!
export boxsmooth!
export outest_bounds #

"""
    outest_bounds(cx,sx) -> px0

Helper function to find maximum padding in pixels required to accomodate all query points `cx` outside of the image size 1:`sx`.

# Arguments:
- `cx`: list of integer star centers (in either x or y)
- `sx`: image dimension along the axis indexed by `cx`

# Outputs:
- `px0`: maximum padding in pixels required to accomodate all query points
"""
function outest_bounds(cx,sx)
    px0 = 0
    sortcx = sort(cx)
    if sortcx[1] < 1
        px0 = abs(sortcx[1]-1)
    end
    if sortcx[end] > sx
        if px0 < (sortcx[end]-sx)
            px0 = (sortcx[end]-sx)
        end
    end
    return px0
end

"""
    boxsmooth!(out::AbstractArray, arr::AbstractArray, tot::Array{T,1}, widx::Int, widy::Int)

Boxcar smooths an input image (or paddedview) `arr` with window size `widx` by
`widy`. We pass the original image size `sx` and `sy` to help handle image views.

# Arguments:
- `out::AbstractArray`: preallocated output array for the boxcar smoothed image
- `arr::AbstractArray`: input array for which boxcar smoothing is computed (generally paddedview)
- `tot::Array{T,1}`: preallocated array to hold moving sums along 1 dimension
- `widx::Int`: size of boxcar smoothing window in x
- `widy::Int`: size of boxcar smoothing window in y
"""
function boxsmooth!(out::AbstractArray, arr::AbstractArray, tot::Array{T,1}, widx::Int, widy::Int) where T
    (sx, sy) = size(arr)
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    for j=1:(sy-widy+1)
        if (j==1)
            for n = 1:widy
                @simd for m = 1:sx
                    @inbounds tot[m] += arr[m,n]
                end
            end
        else
            @simd for m = 1:sx
                @inbounds tot[m] += arr[m,j+widy-1]-arr[m,j-1]
            end
        end
        tt=zero(eltype(out))
        for i=1:(sx-widx+1)
            if (i==1)
                @simd for n=1:widx
                    @inbounds tt += tot[n]
                end
            else
                @inbounds tt += tot[i+widx-1]-tot[i-1]
            end
            @inbounds out[i,j] = tt
        end
    end
end

"""
    cov_avg!(bimage, ism, bism, in_image; Np::Int=33, widx::Int=129, widy::Int=129, ftype::Int=32)

Key function for constructing the (shifted and multiplied) versions of the input image used to quickly
estimate the local covariance matrix at a large number of locations. The main output is in the preallocated
`bism` which is used as an input to `build_cov!`.

# Arguments:
- `bimage`: preallocated output array for the boxcar smoothed unshifted image
- `ism`: preallocated intermediate array for the input image times itself shifted
- `bism`: preallocated output array to store boxcar-smoothed image products for all shifts
- `in_image`: input image the local covariance of which we want to estimate

# Keywords:
- `Np::Int`: size of local covariance matrix in pixels (default 33)
- `widx::Int`: width of boxcar window in x which determines size of region used for samples for the local covariance estimate (default 129)
- `widy::Int`: width of boxcar window in y which determines size of region used for samples for the local covariance estimate (default 129)
- `ftype::Int`: determine the Float precision, 32 is Float32, otherwise Float64
"""
function cov_avg!(bimage, ism, bism, in_image; Np::Int=33, widx::Int=129, widy::Int=129, ftype::Int=32)
    if ftype == 32
        T = Float32
    else
        T = Float64
    end

    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    halfNp = (Np-1) ÷ 2

    (sx1, sy1) = size(in_image)
    tot = zeros(T,sx1);
    boxsmooth!(bimage,in_image,tot,widx,widy)
    # loop over shifts
    for dc=0:Np-1       # column shift loop
        for dr=1-Np:Np-1   # row loop, incl negatives
            if (dr < 0) & (dc == 0)
                continue
            end
            # ism = image, shifted and multipled
            @inbounds ism .= in_image .* ShiftedArrays.circshift(in_image,(-dr, -dc))
            fill!(tot,0)
            boxsmooth!(view(bism,:,:,dr+Np,dc+1),ism,tot,widx,widy) # bism = boxcar(ism)
        end
    end
    return
end
