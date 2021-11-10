using LinearAlgebra

export stamp_cutter
export gen_pix_mask
export condCovEst_wdiag
export build_cov!

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
function stamp_cutter(cx,cy,residimIn,star_im,maskim;Np=33)
    radNp = (Np-1)÷2

    cov_stamp = cx-radNp:cx+radNp,cy-radNp:cy+radNp
    @views data_in = residimIn[cov_stamp[1],cov_stamp[2]]
    @views stars_in = star_im[cov_stamp[1],cov_stamp[2]]
    @views kmasked2d = maskim[cov_stamp[1],cov_stamp[2]];
    return data_in, stars_in, kmasked2d
end

function gen_pix_mask(kmasked2d,psfmodel,circmask,x_star,y_star,flux_star;Np=33,thr=20)

    psft = psfmodel(x_star,y_star,Np)

    if flux_star < 1e4  #these are the pixels we want the cov of
        kpsf2d = (psft .> thr/1e4)
    else
        kpsf2d = (psft .> thr/flux_star)
    end

    kstar = (kmasked2d .| kpsf2d .| circmask)[:]
    cntks = count(kstar)

    dnt = 0
    if cntks > 33^2-248
        dnt = 1
        kstar = (kmasked2d .| kpsf2d)[:]
    end

    if cntks > 33^2-248 #this is about a 10% cut, and is the sum of bndry
        dnt = 2
        kmasked2d[1:2,:] .= 0
        kmasked2d[(end-1):end,:] .= 0
        kmasked2d[:,1:2] .= 0
        kmasked2d[:,(end-1):end] .= 0

        kpsf2d[1:2,:] .= 0
        kpsf2d[(end-1):end,:] .= 0
        kpsf2d[:,1:2] .= 0
        kpsf2d[:,(end-1):end] .= 0

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
function condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=33,export_mean=false,n_draw=0)
    k = .!km
    kstar = kpsf2d[:]
    for i=1:Np*Np cov_loc[i,i] += stars_in[i] end
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

    @views resid_mean = (uncond_input[kstar]'*ipcovCp)./var_wdb
    @views pred_mean = (kstarpred'*ipcovCp)./var_wdb

    # Currently limited to the Np region. Often useful to have some context with a larger
    # surrounding region... TO DO to implement
    out = []
    push!(out,[sqrt(var_wdb^(-1)) resid_mean+pred_mean resid_mean pred_mean])
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
                t = bism[pr+Δr,pc+Δc,dr+Np,dc+1]/(widx*widy) - μ1μ2
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
