export fil_array!
export add_rndm_stars!
export run_fil

using LinearAlgebra
using StatsBase
using cloudCovErr
using PyCall
import Conda

function __init__()
    if !haskey(Conda._installed_packages_dict(),"crowdsourcephoto")
        Conda.add("nomkl",channel="conda-forge")
        Conda.add("crowdsourcephoto",channel="conda-forge")
    end
    py"""
    import crowdsource
    """
end

function fil_array!(image;msz=59)
    (Nx, Ny) = size(image)
    x=1:Nx
    xx = ones(Ny)' .* x
    x0 = 520-msz
    for lam0 in [4,20,50]
        x0 += maximum([6*lam0,2*msz])
        image += 30 .*exp.(-(xx.-x0).^2 ./(lam0^2))

        x0 += maximum([6*lam0,2*msz])
        image -= 30 .*exp.(-(xx.-x0).^2 ./(lam0^2))
    end
    return image
end

function add_rndm_stars!(image,psf_sp,psfmask_sp,Nstar;msz = 59)
    mask = zeros(Bool,size(image))
    cxlst = zeros(Nstar)
    cylst = zeros(Nstar)

    for i=1:Nstar
        rng=MersenneTwister(10+i)
        cx = rand(rng,1:2046-msz)
        cy = rand(rng,1:4094-msz)
        cxlst[i] = cx
        cylst[i] = cy
        mask[cx+1:(cx+msz),cy+1:(cy+msz)] .+= psfmask_sp
        image[cx+1:(cx+msz),cy+1:(cy+msz)] .+=  30000 .*psf_sp
    end
    return mask
end

function run_fil(base,date,filt,vers,basecat,ccd;corrects7=true,msz=59,ftype::Int=32)
    hmsz = (msz-1)÷2

    if ftype == 32
        T = Float32
    else
        T = Float64
    end

    # loads from disk
    ref_im, d_im = read_decam(base,date,filt,vers,ccd,corrects7=corrects7)
    d_im = nothing
    (sx0, sy0) = size(ref_im)
    nodeblend_maskbit = 2^30
    sharp_maskbit = 2^31
    d_im = zeros(Int32,sx0, sy0) .+ nodeblend_maskbit .+ sharp_maskbit;
    med = StatsBase.median(ref_im)
    ref_im = nothing

    in_img = med.*ones(T, sx0, sy0)
    test_out = fil_array!(in_img);

    psfmodel = cloudCovErr.load_psfmodel_cs(basecat,date,filt,vers,ccd)
    psfstatic59 = psfmodel(sx0÷2,sy0÷2,msz)
    psf_mask = zeros(Bool,msz,msz)
    psf_mask[1+hmsz-4:1+hmsz+4,1+hmsz-4:1+hmsz+4] .= 1;

    mask_out = add_rndm_stars!(test_out,psfstatic59,psf_mask,4000);
    add_noise!(test_out,gain;seed=2021)

    in_wts = gain./test_out;

    testres = py"crowdsource.fit_im"(test_out,psfmodel,dq=d_im,weight=sqrt.(in_wts),
        ntilex=4, ntiley=2, refit_psf=false,psfderiv=true,maxstars=320000,verbose=true,threshold=200)

    temp = zeros(size(testres[1])[1],20)
    for i=1:size(testres[1])[1]
        temp[i,:] .= hcat(testres[1][i]...)[:]
    end
    mod_im = testres[2]
    sky_im = testres[3]
    psf_out = testres[4];

    x_stars = temp[:,1]
    y_stars = temp[:,2]
    flux_stars = temp[:,3]
    dflux_stars = temp[:,7]

    return x_stars, y_stars, flux_stars, dflux_stars, mod_im, sky_im, psf_out
end
