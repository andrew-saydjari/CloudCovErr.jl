## Handler for reading outputs of crowdsource processing on DECaPS








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
