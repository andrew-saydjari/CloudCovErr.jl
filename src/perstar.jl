using LinearAlgebra

export stamp_cutter
export gen_pix_mask
export condCovEst_wdiag
export build_cov!

"""
    stamp_cutter(cx,cy,residimIn,star_im,maskim;Np=33) -> data_in, stars_in, kmasked2d

Cuts out local stamps around each star of the various input images to be used for
per star statistics calculations.

# Arguments:
- `cx`: center coorindate x of the stamp
- `cy`: center coorindate y of the stamp
- `residimIn`: residual image with infilling from which covariance was estimated
- `star_im`: input image of model of stars only. abs(mod_im-sky_im)
- `maskim`: input image of upstream masked pixels

# Keywords:
- `Np`: size of covariance matrix footprint around each star (default 33)

# Output:
- `data_in`: local stamp of the (non-infilled) residual image
- `stars_in`: local stamp of model of stars only
- `kmasked2d`: local stamp of upstream masked pixels
"""
function stamp_cutter(cx,cy,residimIn,star_im,maskim;Np=33)
    radNp = (Np-1)÷2

    cov_stamp = cx-radNp:cx+radNp,cy-radNp:cy+radNp
    @views data_in = residimIn[cov_stamp[1],cov_stamp[2]]
    @views stars_in = star_im[cov_stamp[1],cov_stamp[2]]
    @views kmasked2d = maskim[cov_stamp[1],cov_stamp[2]];
    return data_in, stars_in, kmasked2d
end

"""
    gen_pix_mask(kmasked2d,psfmodel,circmask,x_star,y_star,flux_star;Np=33,thr=20) -> psft, kstar[:], kpsf2d, kcond0, kcond, kpred, dnt

Assigns pixels in the local subimage around a star to either be "good", "hidden", or "ignored"
based on user settings and the flux of the star. Reads in masked pixels from the quality
flags on pixels coming from the community pipeline `kmasked2d`, a PSF model for the star,
and a precomputed circular mask `circmask` to exclude pixels at a large radius from the
stellar center since they have little impact on the regression of hidden pixels. The pixels
assigned as "hidden" and to be interpolated are determined by a `thr` on the pixel values
for `flux_star` times the PSF model. We use a parametric PSFs that varies with position and
query the PSF at the stellar position for each star.

# Arguments:
- `kmasked2d`: Bool mask from upstream pixel quality flags to assign pixels as "ignored"
- `psfmodel`: parametric PSF model that can be queried at different positions
- `circmask`: static Bool mask assigning pixels beyond some radius of the stellar center as "ignored"
- `x_star`: x-coordinate of the star (used only for flexible PSF model query)
- `y_star`: y-coordinate of the star (used only for flexible PSF model query)
- `flux_star`: flux of star in ADU to determine how large a region to make "hidden"

# Keywords:
- `Np`: size of local covariance matrix in pixels (default 33)
- `thr`: threshold for psf-based masking of the residuals (larger more "hidden")

# Output:
- `psft`: static array (image) of the stellar PSF
- `kstar`: Boolean indexes the NOT "good" pixels
- `kpsf2d`: Boolean indexes the "hidden" pixels
- `kcond0`: initial number of "good" pixels
- `kcond`: final number of "good" pixels after fallbacks
- `kpred`: the number of pixels "hidden"
- `dnt`: quality flag bits on the solution
"""
function gen_pix_mask(kmasked2d,psfmodel,circmask,x_star,y_star,flux_star;Np=33,thr=20)

    psft = psfmodel(x_star,y_star,Np)

    if flux_star < 1e4  #these are the pixels we want the cov of
        kpsf2d = (psft .> thr/1e4)
    else
        kpsf2d = (psft .> thr/flux_star)
    end

    kstar = (kmasked2d .| kpsf2d .| circmask)
    kcond0 = Np^2-count(kstar)

    dnt = 0
    if kcond0 < 4*Np
        dnt += 1
        kstar = (kmasked2d .| kpsf2d)
    end

    kcond = Np^2-count(kstar)

    if kcond < 4*Np #this is about a 10% cut, and is the sum of bndry
        dnt += 2
        kstar[1,:] .= 0
        kstar[end,:] .= 0
        kstar[:,1] .= 0
        kstar[:,end] .= 0

        kpsf2d[1,:] .= 0
        kpsf2d[end,:] .= 0
        kpsf2d[:,1] .= 0
        kpsf2d[:,end] .= 0
    end

    psfr = minimum(psft)/maximum(psft)
    if psfr < 0
        dnt += 8
        psft .= abs.(psft)
    end
    if psfr < -1e-3
        dnt += 16
    end
    if psfr < -1e-1
        dnt += 32
    end

    kpred = count(kpsf2d)
    kcond = Np^2-count(kstar)

    return psft, kstar[:], kpsf2d, kcond0, kcond, kpred, dnt
end

"""
    condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=33,export_mean=false,n_draw=0,diag_on=true) -> out

Using a local covariance matrix estimate `cov_loc` and a set of known ("good") pixels `km`
and "hidden" pixels `kpsf2d`, this function computes a prediction for the mean value
of the `kpsf2d` pixels and the covariance matrix of the `kpsf2d` pixels. In terms of
statistics use to adjust the photometry of a star, we are only interested in the
pixels masked as a result of the star (i.e. not a detector defect or cosmic ray nearby).
The residual image `data_in` and a model of the counts above the background coming from the
star `stars_in` for the local patch are also inputs of the function. Correction factors for
the photometric flux and flux uncertainities are outputs as well as a chi2 value for the
"good" pixels. The output list can conditionally include the mean reconstruction and
draws from the distribution of reconstructions.

# Arguments:
- `cov_loc`: local covariance matrix
- `μ`: vector containing mean value for each pixel in the patch
- `km`: unmasked pixels
- `kpsf2d`: pixels masked due to the star of interest
- `data_in`: (non-infilled) residual image in local patch
- `psft`: static array (image) of the stellar PSF

# Keywords:
- `Np`: size of local covariance matrix in pixels (default 33)
- `export_mean`: when true, returns the mean conditional prediction for the "hidden" pixels (default false)
- `n_draw`: when nonzero, returns that number of realizations of the conditional infilling (default 0)
- `diag_on`: flag for adding to the pixelwise uncertainty based on the photoelectron counts of the modeled star (default true)

# Output:
- `out[1]`: flux uncertainity of the star
- `out[2]`: flux uncertainity of the star assuming the covariance matrix were diagonal
- `out[3]`: flux correction which must be added to correct the input flux estimate
- `out[4]`: flux correction coming from the residuals (fdb_res)
- `out[5]`: flux correction coming from the predicted background (fdb_pred)
- `out[6]`: chi2 for the "good" pixels under `cov_loc` as a metric on how good our assumptions are
- `out[7]`: local region (image) with "hidden" pixels replaced by the mean conditional estimate (optional output)
- `out[end:end+n_draw]`: local region (image) with "hidden" pixels replaced by the draws from the conditional distribution (optional output)
"""
function condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=33,export_mean=false,n_draw=0,diag_on=true)
    k = .!km
    kstar = kpsf2d[:]
    if diag_on
        for i=1:Np*Np cov_loc[i,i] += stars_in[i] end
    end
    cov_kk = Symmetric(cov_loc[k,k])
    cov_kkstar = cov_loc[k,kstar];
    cov_kstarkstar = cov_loc[kstar,kstar];
    icov_kkC = cholesky(cov_kk)
    icovkkCcovkkstar = icov_kkC\cov_kkstar
    predcovar = Symmetric(cov_kstarkstar - (cov_kkstar'*icovkkCcovkkstar))
    ipcovC = cholesky(predcovar)

    @views uncond_input = data_in[:]
    @views cond_input = data_in[:].- μ

    kstarpredn = (cond_input[k]'*icovkkCcovkkstar)'
    kstarpred = kstarpredn .+ μ[kstar]
    @views p = psft[kpsf2d][:]
    ipcovCp = ipcovC\p

    #@views std_wdiag = sqrt(abs(sum((pw.^(2)).*diag(predcovar[kpsf1d_kstar,kpsf1d_kstar]))))/sum(p2w)
    @views var_wdb = (p'*ipcovCp)

    var_diag = 0
    for i in 1:size(predcovar)[1]
        var_diag += (p[i]^2)/predcovar[i,i]
    end

    @views resid_mean = (uncond_input[kstar]'*ipcovCp)./var_wdb
    @views pred_mean = (kstarpred'*ipcovCp)./var_wdb

    #if we can afford it, a nice check would be how good of a covariance matrix
    #cov is for the k pixels (which we think are clean)
    chi2 = cond_input[k]'*(icov_kkC\cond_input[k])

    # Currently limited to the Np region. Often useful to have some context with a larger
    # surrounding region... TO DO to implement
    out = []
    push!(out,[sqrt(var_wdb^(-1)) sqrt(var_diag^(-1)) pred_mean-resid_mean resid_mean pred_mean chi2])
    if export_mean
        mean_out = copy(data_in)
        mean_out[kstar] .= kstarpred
        push!(out,mean_out)
    end
    if n_draw != 0
        covsvd = svd(predcovar)
        sqrt_cov = covsvd.V*diagm(sqrt.(covsvd.S))*covsvd.Vt;
        noise = sqrt_cov*randn(size(sqrt_cov)[1],n_draw)

        draw_out = repeat(copy(data_in)[:],outer=[1 n_draw])
        draw_out[kstar,:] .= repeat(kstarpred,outer=[1 n_draw]) .+ noise
        push!(out,draw_out)
    end

    return out
end

"""
    build_cov!(cov::Array{T,2},μ::Array{T,1},cx::Int,cy::Int,bimage::Array{T,2},bism::Array{T,4},Np::Int,widx::Int,widy::Int) where T <:Union{Float32,Float64}

Constructs the local covariance matrix and mean for an image patch of size `Np` x `Np` pixels around a location
of interest (`cx`,`cy`). The construction is just a lookup of pixel values from the stored boxcar-smoothed copies
of the input image times itself shifted in `bism`. Passing the smoothed image `bimage` and the widths of the boxcar
mean `widx` and `widy` is helpful for the mean and normalization. The covariance and mean are updated in place
for speed since this operation may be performed billions of times since we construct a new covariance matrix for
every detection. Math may either be performed `Float32` or `Float64`.

# Arguments:
- `cov::Array{T,2}`: preallocated output array for local covariance matrix
- `μ::Array{T,1}`: preallocated output vector for local mean
- `cx::Int`: x-coordinate of the center of the local region
- `cy::Int`: y-coordinate of the center of the local region
- `bimage::Array{T,2}`: boxcar smoothed unshifted image
- `bism::Array{T,4}`: boxcar-smoothed image products for all shifts
- `Np::Int`: size of local covariance matrix in pixels
- `widx::Int`: width of boxcar window in x which determines size of region used for samples for the local covariance estimate
- `widy::Int`: width of boxcar window in y which determines size of region used for samples for the local covariance estimate
"""
function build_cov!(cov::Array{T,2},μ::Array{T,1},cx::Int,cy::Int,bimage::Array{T,2},bism::Array{T,4},Np::Int,widx::Int,widy::Int) where T <:Union{Float32,Float64}
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    halfNp = (Np-1) ÷ 2
    Δr, Δc = cx-(halfNp+1), cy-(halfNp+1)
    # Δr, Δc = cx-(halfNp-1), cy-(halfNp-1)
    for dc=0:Np-1       # column shift loop
        pcr = 1:Np-dc
        for dr=1-Np:Np-1# row loop, incl negatives
            if (dr < 0) & (dc == 0)
                continue
            end
            if dr >= 0
                prr = 1:Np-dr
            end
            if (dr < 0) & (dc > 0)
                prr = 1-dr:Np
            end

            for pc=pcr, pr=prr
                i = ((pc   -1)*Np)+pr
                j = ((pc+dc-1)*Np)+pr+dr
                @inbounds μ1μ2 = bimage[pr+Δr,pc+Δc]*bimage[pr+dr+Δr,pc+dc+Δc]/((widx*widy)^2)
                @inbounds t = bism[pr+Δr,pc+Δc,dr+Np,dc+1]/(widx*widy) - μ1μ2
                @inbounds cov[i,j] = t
                @inbounds cov[j,i] = t
                if i == j
                    @inbounds μ[i] = μ1μ2
                end
            end
        end
    end
    cov .*= (widx*widy)/((widx*widy)-1)
    return
end
