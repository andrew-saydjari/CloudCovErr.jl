# API Reference

## Covariance Construction Functions

```@docs
cov_avg!
boxsmooth!
outest_bounds
```

## Per Star Functions

```@docs
stamp_cutter
gen_pix_mask
condCovEst_wdiag
build_cov!
```

## Image Infill and Masking

```@docs
prelim_infill!
gen_mask_staticPSF!
gen_mask_staticPSF2!
im_subrng
add_sky_noise!
add_noise!
findmaxpsf
kstar_circle_mask
```

## DECam Specific Functions

```@docs
read_decam
read_crowdsource
inject_rename
load_psfmodel_cs
save_fxn
get_catnames
proc_ccd
proc_all
```
