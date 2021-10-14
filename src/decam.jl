## Handler for reading outputs of crowdsource processing on DECaPS

#

date = "170420_040428"
key = "N14"
filt = "g"

base = "/n/fink2/decaps/c4d_"
vers = "v1"

f = FITS(base*date*"_ooi_"*filt*"_"*vers*".fits.fz")
ref_im = read(f[key])
close(f)
f = FITS(base*date*"_ood_"*filt*"_"*vers*".fits.fz")
d_im = read(f[key])
close(f)

f = FITS(base*date*"_oow_"*filt*"_"*vers*".fits.fz")
w_im = read(f[key])
close(f)

f = FITS("/n/home12/saydjari/finksagescratch/decaps/cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits")
x_stars = read(f[key*"_CAT"],"x")
y_stars = read(f[key*"_CAT"],"y")
flux_stars = read(f[key*"_CAT"],"flux")
dflux_stars = read(f[key*"_CAT"],"dflux")
fracflux_stars = read(f[key*"_CAT"],"fracflux")
flags_stars = read(f[key*"_CAT"],"flags")
pN_stars = read(f[key*"_CAT"],"prN")
pL_stars = read(f[key*"_CAT"],"prL")
pR_stars = read(f[key*"_CAT"],"prR")
pE_stars = read(f[key*"_CAT"],"prE")
pno_stars = read(f[key*"_CAT"],"passno")
close(f)

f = FITS("/n/home12/saydjari/finksagescratch/decaps/mod/c4d_"*date*"_ooi_"*filt*"_"*vers*".mod.fits")
mod_im = read(f[key*"_MOD"])
sky_im = read(f[key*"_SKY"])
msk_im = read(f[key*"_MSK"])
close(f)

diffuse = flags_stars .& 2^21
brightstar = flags_stars .& 2^23
galaxy = flags_stars .& 2^24

# Plotting WU
fig = plt.figure(figsize=(16,8), dpi=150)
plt.subplots_adjust(wspace=0.2,hspace=0.1)
plt.suptitle(date*" "*filt*" "*key,y=0.91,fontsize=14)
ax = fig.add_subplot(2,2,1)
input = ref_im .-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.scatter(x_stars,y_stars,s=1,color="red", edgecolor="none",alpha=1)
ax.set_title("Image")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

ax = fig.add_subplot(2,2,2)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.scatter(x_stars,y_stars,s=1,c=diffuse,edgecolor="none",alpha=1,cmap="cet_bkr_r")
ax.set_title("Nebulous Mask")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

ax = fig.add_subplot(2,2,3)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.scatter(x_stars,y_stars,s=1,c=galaxy,edgecolor="none",alpha=1,cmap="cet_bkr_r")
ax.set_title("Galaxy Mask")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

ax = fig.add_subplot(2,2,4)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.scatter(x_stars,y_stars,s=1,c=brightstar,edgecolor="none",alpha=1,cmap="cet_bkr_r")
ax.set_title("Bright Stars Mask")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

#
fig = plt.figure(figsize=(16,8), dpi=150)
plt.subplots_adjust(wspace=0.2,hspace=0.1)
plt.suptitle(date*" "*filt*" "*key,y=0.91,fontsize=14)
ax = fig.add_subplot(2,2,1)
input = ref_im .-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
sc = ax.scatter(x_stars,y_stars,s=2,c=pN_stars,edgecolor="none",alpha=1,cmap="cet_bgy")
ax.set_title("Nebulosity Probability")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc, cax=cax)

ax = fig.add_subplot(2,2,2)
input = ref_im .-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
sc = ax.scatter(x_stars,y_stars,s=2,c=pL_stars,edgecolor="none",alpha=1,cmap="cet_bgy")
ax.set_title("Nebulosity Light Probability")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc, cax=cax)

ax = fig.add_subplot(2,2,3)
input = ref_im .-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
sc = ax.scatter(x_stars,y_stars,s=2,c=pR_stars,edgecolor="none",alpha=1,cmap="cet_bgy")
ax.set_title("Regular Probability")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc, cax=cax)

ax = fig.add_subplot(2,2,4)
input = ref_im .-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
sc = ax.scatter(x_stars,y_stars,s=2,c=pE_stars,edgecolor="none",alpha=1,cmap="cet_bgy")
ax.set_title("Error Probability")
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc, cax=cax)

#
x_limit = 700
y_limit = 600
sizef = 300
sp = 2

fig = plt.figure(figsize=(16,8), dpi=150)
fig.subplots_adjust(wspace=0.1, hspace=0.1)
plt.suptitle("Subimage "*date*" "*filt*" "*key,y=0.965,fontsize=14)
#First the reference info
ax = fig.add_subplot(2,4,5)
input = ref_im.-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Raw Image")

ax = fig.add_subplot(2,4,2)
input = msk_im
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1",
    aspect="equal",
    vmin=0,
    vmax=1
)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Neb Mask")


ax = fig.add_subplot(2,4,3)
ax.imshow(
    w_im,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Weights")

ax = fig.add_subplot(2,4,4)
ax.imshow(
    d_im,
    origin="lower",
    interpolation="nearest",
    cmap="cet_fire",
    aspect="equal",
)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Bitmask")

ax = fig.add_subplot(2,4,1)
input = ref_im.-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.plot(x_stars,y_stars,"r.",ms=sp)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Stars")

ax = fig.add_subplot(2,4,6)
input = mod_im.-median(mod_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)

ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Model")

ax = fig.add_subplot(2,4,7)
resid = mod_im.-ref_im
input = resid

ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_bkr",
    aspect="equal",
    vmin = -50,
    vmax = 50
)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Residuals")

ax = fig.add_subplot(2,4,8)
input = sky_im.-median(sky_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Background");

#
sizef = 500

fig = plt.figure(figsize=(12,8), dpi=150)
fig.subplots_adjust(wspace=0.2, hspace=0.1)
plt.suptitle(date*" "*filt*" "*key, y=0.82,fontsize=14)

ax = fig.add_subplot(1,2,1)
maxiter = maximum(pno_stars)
cmap = plt.get_cmap("tab10", Int(maxiter+1))
input = ref_im.-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
sc = ax.scatter(x_stars,y_stars,s=5,c=pno_stars,cmap=cmap, edgecolor="none",alpha=1)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("Deblending Order")

norm2 = plt.matplotlib.colors.Normalize(vmin=0,vmax=maxiter+1)
cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
cbar = ax.figure.colorbar(
            plt.matplotlib.cm.ScalarMappable(norm=norm2, cmap=cmap),
            #boundaries=0:10-0.5,
            ticks=collect(-0.5 .+ (0:maxiter+1)),
            cax=cax)#
cbar.ax.set_yticklabels(collect(-1:maxiter));

ax = fig.add_subplot(1,2,2)
input = ref_im.-median(ref_im)
ax.imshow(
    input,
    origin="lower",
    interpolation="nearest",
    cmap="cet_CET_L1_r",
    aspect="equal",
    vmin=-50,
    vmax=50
)
sc = ax.scatter(x_stars,y_stars,s=5,c=fracflux_stars,cmap="cet_bgy", vmin=0,vmax=1, edgecolor="none",alpha=1)
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.set_xlim((x_limit,x_limit+sizef))
ax.set_ylim((y_limit,y_limit+sizef))
ax.set_title("FracFlux")

cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc, cax=cax)

#
f = FITS("/n/home12/saydjari/finksage/Working/2021_10_07/cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits")
x_stars0 = read(f[end-3],"x")
y_stars0 = read(f[end-3],"y")
flux_stars0 = read(f[end-3],"flux")
x_stars = read(f[end],"x")
y_stars = read(f[end],"y")
flux_stars = read(f[end],"flux")
dflux_stars = read(f[end],"dflux")
key = split(read_header(f[end])["EXTNAME"],"_")[1]
close(f)

#
i=1
l=0
d=0

merged_cat = []
merged_ideal_cat = []
merged_ideal_cat_full = []
detected_mask = []

d_cutoff = 3^2

key_ind=0
for k=1:size(x_stars0)[1]
    dist = (x_stars.-x_stars0[k]).^2 .+ (y_stars.-y_stars0[k]).^2
    ind = dist .< d_cutoff
    if count(ind) == 1
        pos = findall(ind)
        push!(detected_mask,1)
        push!(merged_cat, vec([key_ind i k x_stars[ind] y_stars[ind] flux_stars[ind] dflux_stars[ind] dist[ind] pos]))
        push!(merged_ideal_cat, vec([key_ind i k x_stars0[k] y_stars0[k] flux_stars0[k]]))
        push!(merged_ideal_cat_full, vec([key_ind i k x_stars0[k] y_stars0[k] flux_stars0[k]]))
    elseif count(ind) > 1
        #println("2 matches")
        push!(detected_mask,0)
        push!(merged_ideal_cat_full, vec([key_ind i k x_stars0[k] y_stars0[k] flux_stars0[k]]))
        d+=1
    else
        #println("no matches")
        push!(detected_mask,0)
        push!(merged_ideal_cat_full, vec([key_ind i k x_stars0[k] y_stars0[k] flux_stars0[k]]))
        l+=1
    end
    i+=1
end

detected_mask1 = convert(BitArray,detected_mask)
merged_ideal_cat_full1 = hcat(merged_ideal_cat_full...)
merged_ideal_cat1 = hcat(merged_ideal_cat...)
merged_cat1 = hcat(merged_cat...)

println("Total Num Stars: "*string(i-1))
println("No Matches: "*string(l))
println("Too Many Matches: "*string(d))
println("Total Good: "*string(i-1-d-l))

#
fig = plt.figure(figsize=(24,16), dpi=150)
fig.subplots_adjust(wspace=0.3, hspace=0.2)
lperc = 2
hperc = 97

ax = fig.add_subplot(2,3,1)
plt.suptitle(date*" "*filt*" "*key, y=0.92,fontsize=14)
maskpos = flux_stars.>0;

data = (-2.5*log10.(flux_stars[maskpos]).+zpdic[filt])
ax.hist(data,bins=50,alpha=0.5,range=(16.5,24.5));
ax.set_xlabel("Mags")
ax.set_title("Distribution of Stellar Flux")

ax = fig.add_subplot(2,3,4)
maskpos = merged_cat1[6,:].>0;

data = (-2.5*log10.(merged_cat1[6,maskpos]).+zpdic[filt])
ax.hist(data,bins=50,alpha=0.5,range=(16.5,24.5),histtype="step");

maskpos1 = merged_ideal_cat1[6,:].>0;

data = (-2.5*log10.(merged_ideal_cat1[6,maskpos1]).+zpdic[filt])
ax.hist(data,bins=50, alpha=0.5, range=(16.5,24.5),histtype="step");

ax.legend(["Recovered", "Inserted"],loc="upper left")
ax.set_xlabel("Mags")

diffy = (merged_cat1[6,:].-merged_ideal_cat1[6,:]);
zdiff = 100 .*diffy./maximum(hcat(merged_cat1[6,:],merged_ideal_cat1[6,:]),dims=2);
fluxiqr = StatsBase.iqr(zdiff[:])
lower, upper = (-10,10)

maskpos = (merged_cat1[6,:].>0) .& (merged_ideal_cat1[6,:].>0);

ax = fig.add_subplot(2,3,2)
ax.hist(zdiff[:],bins=100,alpha=0.5,density=true,range=(lower, upper));
ax.set_xlim([lower,upper])
ax.set_title("Distribution of Percent Errors")
ax.set_xlabel("Percent Diff")
ax.set_ylabel("Density")
ax.legend([" IQR = "*string(round(fluxiqr,digits = 2))*" \n Std = "*string(round(StatsBase.std(zdiff),digits = 2))*" \n Med = "*string(round(median(zdiff),digits = 2))])

ax = fig.add_subplot(2,3,5)
sc = ax.hist2d(zdiff[maskpos], -2.5*log10.(merged_cat1[6,maskpos]).+zpdic[filt],
    bins=[100,100],
    range=[(lower,upper),(16.5,24.5)],
    cmap="cet_fire",
)
ax.set_ylabel("Observed Mags")
ax.set_xlabel("Percent Difference")
cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc[4], cax=cax)

diffy_a = diffy .- median(zdiff)/100 .*maximum(hcat(merged_cat1[6,:],merged_ideal_cat1[6,:]),dims=2)
errbars_floor = sqrt.(merged_cat1[7,:].^2 .+ 0.0001.*merged_cat1[6,:].^2)
zdiffo = diffy_a./errbars_floor
nan_mask_o = .!isnan.(zdiffo)
lowerp, upperp = StatsBase.percentile(zdiffo[nan_mask_o],[lperc,hperc])

ax = fig.add_subplot(2,3,3)
ax.hist(zdiffo[nan_mask_o],bins=100,alpha=0.5,density=true,range=(lowerp, upperp));
ax.set_xlim([lowerp,upperp])
ax.set_title("Distribution of Z-Scores")
ax.set_xlabel("Z-Score")
ax.set_ylabel("Density")
ax.legend([" IQR = "*string(round(StatsBase.iqr(zdiffo[nan_mask_o]),digits = 2))*" \n Std = "*string(round(StatsBase.std(zdiffo[nan_mask_o]),digits = 2))*" \n Med = "*string(round(median(zdiffo[nan_mask_o]),digits = 2))])

ax = fig.add_subplot(2,3,6)
sc = ax.hist2d(zdiffo[nan_mask_o], -2.5*log10.(merged_cat1[6,maskpos]).+zpdic[filt],
    bins=[100,100],
    range=[(lowerp,upperp),(16.5,24.5)],
    cmap="cet_fire",
)
ax.set_ylabel("Observed Mags")
ax.set_xlabel("Z-Score")
cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(sc[4], cax=cax)

#
thr = 20
data_w = deepcopy(w_im)
rad=255
maskd = (data_w.<=0.0002) .| (d_im .!= 0);
for i=1:size(x_stars)[1]
    fluxt=flux_stars[i]
    if fluxt .> 1e4
        x_star = round(Int64, x_stars[i])
        y_star = round(Int64, y_stars[i])
        mskt = (psf0 .> thr/(fluxt))[maximum([1,1+rad-y_star]):minimum([2046-y_star+rad,511]),maximum([1,1+rad-x_star]):minimum([4094-x_star+rad,511])]
        maskd[maximum([1,1+y_star-rad]):minimum([1+y_star+rad,2046]),maximum([1,1+x_star-rad]):minimum([1+x_star+rad,4094])] .|= mskt
    end
end
resid = deepcopy(mod_im.-ref_im);
resid[maskd] .= NaN;

mask0 =maskd;
mask = convert.(Int64,.!mask0);

testim = deepcopy(mod_im.-ref_im)[550:1400,550:1750];
testim0 = deepcopy(mod_im.-ref_im)[550:1400,550:1750];
skyim = deepcopy(ref_im[550:1400,550:1750].-(mod_im[550:1400,550:1750].-sky_im[550:1400,550:1750]))
skyim0 = deepcopy(ref_im[550:1400,550:1750].-(mod_im[550:1400,550:1750].-sky_im[550:1400,550:1750]))
maskim = deepcopy(mask)[550:1400,550:1750];
maskim0 = deepcopy(mask0)[550:1400,550:1750];
data_w = deepcopy(data_w)[550:1400,550:1750];

skyim[maskim0] .= 0;
testim2 = deepcopy(testim);
skyim2 = deepcopy(skyim);
maskim2 = deepcopy(maskim0);
ref_im0 = ref_im[550:1400,550:1750];
mod_im0 = mod_im[550:1400,550:1750];
gain = median(data_w.*ref_im[550:1400,550:1750])
skyim3 = deepcopy(sky_im)[550:1400,550:1750];

"""
    prelim_infill!(ind)

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

@showprogress for j=1:size(testim2)[2], i=1:size(testim2)[1]
    if maskim0[i,j]
        intermed = -(rand(Distributions.Poisson(convert(Float64,gain*(skyim3[i,j]-testim2[i,j]))))/gain.-skyim3[i,j])
        testim2[i,j] = intermed
    end
end

stars_interior = ((x_stars) .> 550 .+Np) .& ((x_stars) .< 1750 .-Np) .& ((y_stars) .> 550 .+Np) .& ((y_stars) .< 1400 .-Np);
calstar = (x_stars .∈ [merged_cat1[4,:]]) .& stars_interior;
cy = round.(Int64,(x_stars.-549)[calstar]).+1;
cx = round.(Int64,(y_stars.-549)[calstar]).+1;
count(stars_interior), count(calstar)

cov_loc, cnts_loc, μ_loc  = cov_construct_I(testim2,maskim,cy,cx,Np=33,wid=127);

psf33 = py"psf0(0,0,stampsz=33)"
calflux = flux_stars[calstar];

py"""
def load_psfmodel(outfn, key, filter, pixsz=9):
    f = fits.open(outfn)
    psfmodel = psfmod.linear_static_wing_from_record(f[key+"_PSF"].data[0],filter=filter)
    f.close()
    psfmodel.fitfun = partial(psfmod.fit_linear_static_wing, filter=filter, pixsz=pixsz)
    return psfmodel
"""
outfn = "/n/home12/saydjari/finksage/Working/2021_10_07/cat/c4d_"*date*"_ooi_"*filt*"_"*vers*".cat.fits"
key="N14"
psfmodel0 = py"load_psfmodel"(outfn,key,filt)


#
out0 = @showprogress map(per_star,1:size(cx)[1]);

cs_out_0=zeros(size(cx)[1],12)
for i=1:size(cx)[1]
    cs_out_0[i,:]=out0[i][1]
end

stamps1=zeros(size(cx)[1],65,65)
for i=1:size(cx)[1]
    stamps1[i,:,:]=out0[i][2]
end

test2 = merged_cat1[9,:] .∈ [findall(calstar)]
p = sortperm(merged_cat1[9,test2]);

map1 = StatsBase.countmap(merged_cat1[9,test2])
dmask = []
for x in merged_cat1[9,test2][p]
    push!(dmask,map1[x]==1)
end
sum(dmask)

dmask1 = []
for x in(1:7304)[calstar]
    push!(dmask1,map1[x]==1)
end
sum(dmask1)

dmask_1 = convert.(Bool,dmask)
dmask1_1 = convert.(Bool,dmask1);




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
