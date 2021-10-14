## Handler for reading outputs of crowdsource processing on DECaPS
import FITSIO

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
    read_crowdsource(base,date,filt,vers,ccd) -> x_stars, y_stars, decapsid, mod_im, sky_im

Read in outputs of crowdsource, a photometric pipeline. To pair with an arbitrary
photometric pipeline, an analogous read in function should be created. The relevant
outputs are the model image (including the sources) so that we can produce the
residual image, the sky/background model (no sources), and the coordinates of the stars.
The survey id number is also readout of the pipeline solution file to help
cross-validate matching of the disCovErr outputs and the original sources.

# Arguments:
- `base`: parent directory and file name prefix for crowdsource results files
- `date`: date_time of the exposure
- `filt`: optical filter used to take the exposure
- `vers`: NOAO community processing version number
- `ccd`: which ccd we are pulling the image for
"""
function read_crowdsource(base,date,filt,vers,ccd)
    f = FITSIO.FITS("/n/home12/saydjari/finksagescratch/decaps/cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits")
    x_stars = FITSIO.read(f[ccd*"_CAT"],"x")
    y_stars = FITSIO.read(f[ccd*"_CAT"],"y")
    decapsid = FITSIO.read(f[ccd*"_CAT"],"decapsid")
    FITSIO.close(f)

    f = FITSIO.FITS("/n/home12/saydjari/finksagescratch/decaps/mod/c4d_"*date*"_ooi_"*filt*"_"*vers*".mod.fits")
    mod_im = FITSIO.read(f[ccd*"_MOD"])
    sky_im = FITSIO.read(f[ccd*"_SKY"])
    FITSIO.close(f)
    return x_stars, y_stars, decapsid, mod_im, sky_im
end

# #
# thr = 20
# data_w = deepcopy(w_im)
# rad=255
# maskd = (data_w.<=0.0002) .| (d_im .!= 0);
# for i=1:size(x_stars)[1]
#     fluxt=flux_stars[i]
#     if fluxt .> 1e4
#         x_star = round(Int64, x_stars[i])
#         y_star = round(Int64, y_stars[i])
#         mskt = (psf0 .> thr/(fluxt))[maximum([1,1+rad-y_star]):minimum([2046-y_star+rad,511]),maximum([1,1+rad-x_star]):minimum([4094-x_star+rad,511])]
#         maskd[maximum([1,1+y_star-rad]):minimum([1+y_star+rad,2046]),maximum([1,1+x_star-rad]):minimum([1+x_star+rad,4094])] .|= mskt
#     end
# end
# resid = deepcopy(mod_im.-ref_im);
# resid[maskd] .= NaN;
#
# mask0 =maskd;
# mask = convert.(Int64,.!mask0);
#
# testim = deepcopy(mod_im.-ref_im)[550:1400,550:1750];
# testim0 = deepcopy(mod_im.-ref_im)[550:1400,550:1750];
# skyim = deepcopy(ref_im[550:1400,550:1750].-(mod_im[550:1400,550:1750].-sky_im[550:1400,550:1750]))
# skyim0 = deepcopy(ref_im[550:1400,550:1750].-(mod_im[550:1400,550:1750].-sky_im[550:1400,550:1750]))
# maskim = deepcopy(mask)[550:1400,550:1750];
# maskim0 = deepcopy(mask0)[550:1400,550:1750];
# data_w = deepcopy(data_w)[550:1400,550:1750];
#
# skyim[maskim0] .= 0;
# testim2 = deepcopy(testim);
# skyim2 = deepcopy(skyim);
# maskim2 = deepcopy(maskim0);
# ref_im0 = ref_im[550:1400,550:1750];
# mod_im0 = mod_im[550:1400,550:1750];
# gain = median(data_w.*ref_im[550:1400,550:1750])
# skyim3 = deepcopy(sky_im)[550:1400,550:1750];

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
function prelim_infill!(testim,maskim,bimage,bimageI,testim2, maskim2, goodpix; widx = 19, widy=19)

    wid = maximum([widx, widy])
    Δ = (wid-1)÷2
    (sx, sy) = size(testim)

    #the masked entries in testim must be set to 0 so they drop out of the mean
    testim[maskim] .= 0;
    maskim2 .= copy(maskim)
    testim2 .= copy(testim)

    #loop to try masking at larger and larger smoothing to infill large holes
    cnt=0
    while any(maskim2) .& (cnt .< 10)
        in_image = padarray(testim,Pad(:reflect,(Δ+2,Δ+2)));
        in_mask = padarray(.!maskim,Pad(:reflect,(Δ+2,Δ+2)));

        boxsmoothMod!(bimage, in_image, widx, widy, sx, sy)
        boxsmoothMod!(bimageI, in_mask, widx, widy, sx, sy)

        goodpix = (bimageI .> 10)

        testim2[maskim2 .& goodpix] .= (bimage./bimageI)[maskim2 .& goodpix]
        maskim2[goodpix] .= false

        # update loop params
        cnt+=1
        wid*=1.4
        wid = round(Int,wid)
        Δ = (wid-1)÷2
    end
    println((cnt,wid))

    #catastrophic failure fallback
    if cnt == 10
        testim2[maskim2] .= median(in_image)
        println("Infilling Failed Badly")
    end
    return
end

# @showprogress for j=1:size(testim2)[2], i=1:size(testim2)[1]
#     if maskim0[i,j]
#         intermed = -(rand(Distributions.Poisson(convert(Float64,gain*(skyim3[i,j]-testim2[i,j]))))/gain.-skyim3[i,j])
#         testim2[i,j] = intermed
#     end
# end
#
# stars_interior = ((x_stars) .> 550 .+Np) .& ((x_stars) .< 1750 .-Np) .& ((y_stars) .> 550 .+Np) .& ((y_stars) .< 1400 .-Np);
# calstar = (x_stars .∈ [merged_cat1[4,:]]) .& stars_interior;
# cy = round.(Int64,(x_stars.-549)[calstar]).+1;
# cx = round.(Int64,(y_stars.-549)[calstar]).+1;
# count(stars_interior), count(calstar)
#
#### cov_loc, cnts_loc, μ_loc  = cov_construct_I(testim2,maskim,cy,cx,Np=33,wid=127);
#
# psf33 = py"psf0(0,0,stampsz=33)"
# calflux = flux_stars[calstar];
#
# py"""
# def load_psfmodel(outfn, key, filter, pixsz=9):
#     f = fits.open(outfn)
#     psfmodel = psfmod.linear_static_wing_from_record(f[key+"_PSF"].data[0],filter=filter)
#     f.close()
#     psfmodel.fitfun = partial(psfmod.fit_linear_static_wing, filter=filter, pixsz=pixsz)
#     return psfmodel
# """
# outfn = "/n/home12/saydjari/finksage/Working/2021_10_07/cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits"
# key="N14"
# psfmodel0 = py"load_psfmodel"(outfn,key,filt)
#
#
# #
# out0 = @showprogress map(per_star,1:size(cx)[1]);
#
# cs_out_0=zeros(size(cx)[1],12)
# for i=1:size(cx)[1]
#     cs_out_0[i,:]=out0[i][1]
# end
#
# stamps1=zeros(size(cx)[1],65,65)
# for i=1:size(cx)[1]
#     stamps1[i,:,:]=out0[i][2]
# end


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
