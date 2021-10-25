using LinearAlgebra

export stamp_cutter
export gen_pix_mask
export condCovEst_wdiag

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
function condCovEst_wdiag(cov_loc,μ,kstar,kpsf2d,data_in,data_w,stars_in,psft;export_mean=false)
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

    if export_mean
        mean_out = copy(data_in)
        mean_out[kstar] .= kstarpred

        return [std_w std_wdiag sqrt(var_wdb) resid_mean+pred_mean resid_mean pred_mean chi20], mean_out
    end

    return [std_w std_wdiag sqrt(var_wdb) resid_mean+pred_mean resid_mean pred_mean chi20]
end
