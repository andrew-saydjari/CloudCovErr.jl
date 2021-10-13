## utility functions
using OffsetArrays

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
function cov_construct(img, cx, cy; Np::Int=33, widx::Int=129, widy::Int=129)
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    halfNp = (Np-1) ÷ 2
    (Nstar,) = size(cx)
    (sx, sy) = size(img)

    # preallocate output covariance array
    Δr = zeros(Nstar)
    Δc = zeros(Nstar)
    cov = zeros(Nstar,Np*Np,Np*Np)
    μ = zeros(Nstar,Np*Np)

    in_image = padarray(img,Pad(:reflect,(Np+Δx+2,Np+Δy+2)));
    bism = copy(img)
    bimage = copy(img)

    Δr, Δc = cy.-(halfNp-1), cx.-(halfNp-1)

    boxsmoothMod!(in_image,widx,widy,bimage,sx,sy)
    # loop over shifts
    for dc=0:Np-1       # column shift loop
        for dr=1-Np:Np-1   # row loop, incl negatives
            # ism = image, shifted and multipled
            ism = in_image .* OffsetArray(ShiftedArrays.circshift(in_image.parent,(-dr, -dc)), OffsetArrays.Origin(in_image.offsets.+1))

            # bism = boxcar(ism)
            boxsmoothMod!(ism,widx,widy,bism,sx,sy)
            if dr >= 0
                for pc=1:Np-dc, pr=1:Np-dr
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
            if (dr < 0) & (dc > 0)
                ## FIX ME: We can probably remove the boxsmooth compute for dr<0, dc=0 with more care
                for pc=1:Np-dc, pr=1-dr:Np
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
    end
    return (widx*widy)/((widx*widy)-1).*cov, μ # before using, must call Symmetric(cov[i,:,:])
end

"""
    per_star(ind) -> per_star_stats()

    TBD
"""
function per_star(ind)
    Np = 33
    # Mean Handling

    cyy = cy[ind]
    cxx = cx[ind]

    μ = μ_loc[ind,:]
    return per_star_stats(Symmetric(cov_loc[ind,:,:]),cxx,cyy,μ,ind)
end

"""
    per_star_stats(cov_loc,cxx,cyy,μ,ind) -> per_star_stats()

    A bunch o stats
"""
function per_star_stats(cov_loc,cxx,cyy,μ,ind)
    Np = 33
    Npix = 65 #may be able to deprecate the cor stamp during the rewrite
    Nf = 9
    thr = 20

    cy = round(Int64,cyy)
    cx = round(Int64,cxx)

    radNp = (Np-1)÷2
    radNf = (Nf-1)÷2

    Nmp = (Npix-1)÷2
    Nm = Np÷2+1

    cor_stamp = cx-(Nmp):cx+Nmp,cy-(Nmp):cy+Nmp
    cov_stamp = Nmp-radNp+1:Nmp+radNp+1,Nmp-radNp+1:Nmp+radNp+1
    psf_stamp = Nmp-radNf+1:Nmp+radNf+1,Nmp-radNf+1:Nmp+radNf+1
    psf_stamp2 = (radNp)-radNf+1:(radNp)+radNf+1,(radNp)-radNf+1:(radNp)+radNf+1

    @views stars_cut = (mod_im0.-skyim3)[cor_stamp[1],cor_stamp[2]]
    @views stars_in = stars_cut[cov_stamp[1],cov_stamp[2]][:]

    @views mask_cut = maskim0[cor_stamp[1],cor_stamp[2]]
    @views kmasked2d = mask_cut[cov_stamp[1],cov_stamp[2]];
    kmasked2d = mask_cut[cov_stamp[1],cov_stamp[2]] #these are masked from dim or other stars

    psft = psfmodel0(cyy,cxx,33)'

    if calflux[ind] < 1e4  #these are the pixels we want the cov of
        kpsf2d = (psft .> thr/1e4)
    else
        kpsf2d = (psft .> thr/calflux[ind])
    end

    kstar = (kmasked2d .| kpsf2d)[:]
    k = .!kstar;

    cntk = count(k)
    cntks = count(kstar)

    dnt = 0
    if cntk < 128 #this is about a 10% cut, it is actually the sum of the outside bndry in my mind
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
        k = .!kstar;
    end

    #is still the right prefactor for the icov_kk which is the only inverse I really care about?
    ncts = median(Symmetric(cnts_loc[ind,:,:]))
    f = (ncts-sum(k)-2)/(ncts-1)
    f1 = (ncts-length(k)-2)/(ncts-1)

    cov_r = Symmetric(cov_loc) + diagm(stars_in)
    cov_kk = Symmetric(cov_r[k,k])
    cov_kstark = cov_r[kstar,k];
    cov_kstarkstar = Symmetric(cov_r[kstar,kstar]);
    icov_kk = f*inv(cov_kk);
    predcovar = Symmetric(cov_kstarkstar - (cov_kstark*icov_kk*cov_kstark')); #does this need similar icov correct?
    ipcov = inv(predcovar);
    ifull = f1*inv(cov_r);

    @views data_in = testim2[cor_stamp[1],cor_stamp[2]]
    @views uncond_input = data_in[cov_stamp[1],cov_stamp[2]][:]
    @views cond_input = data_in[cov_stamp[1],cov_stamp[2]][:].- μ
    kstarpred = cov_kstark*icov_kk*cond_input[k] .+ μ[kstar]
    kstarpredn = cov_kstark*icov_kk*cond_input[k]

    pdc = isposdef(cov_r)

    chi20 = (kstarpredn'*ipcov*kstarpredn)

    trial = deepcopy(testim[cor_stamp[1],cor_stamp[2]])
    trial[cov_stamp[1],cov_stamp[2]][kpsf2d] .= kstarpred[kpsf2d[:][kstar]];
    @views test_v = trial[cov_stamp[1],cov_stamp[2]][:] .- μ;
    chi22 = (test_v'*ifull*test_v)

    @views data_w_cut = data_w[cov_stamp[1],cov_stamp[2]]
    @views p = psft[kpsf2d][:]
    @views pw = p.*data_w_cut[kpsf2d][:]
    @views p2w = p.*pw

    @views var_w = sqrt(abs(pw'*predcovar[kpsf2d[:][kstar],kpsf2d[:][kstar]]*pw))/sum(p2w)
    @views var_wdiag = sqrt(abs(sum((pw.^(2)).*diag(predcovar[kpsf2d[:][kstar],kpsf2d[:][kstar]]))))/sum(p2w)
    @views var_wdb = (p'*ipcov[kpsf2d[:][kstar],kpsf2d[:][kstar]]*p)

    @views psf_kstar_mean1 = (p'*ipcov[kpsf2d[:][kstar],kpsf2d[:][kstar]]*uncond_input[kpsf2d[:]])./var_wdb

    @views psf_kstar_mean2 = (p'*ipcov[kpsf2d[:][kstar],kpsf2d[:][kstar]]*kstarpred[kpsf2d[:][kstar]])./var_wdb

    @views psf_kstar_mean3 = (p'*ipcov[kpsf2d[:][kstar],kpsf2d[:][kstar]]*kstarpredn[kpsf2d[:][kstar]])./var_wdb

    return [var_w var_wdiag var_wdb psf_kstar_mean1 psf_kstar_mean2 psf_kstar_mean3 chi20 chi22 pdc cntk cntks dnt], trial
end
