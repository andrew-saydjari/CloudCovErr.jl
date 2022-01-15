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
CloudCovErr.decam.read_decam
CloudCovErr.decam.read_crowdsource
CloudCovErr.decam.inject_rename
CloudCovErr.decam.load_psfmodel_cs
CloudCovErr.decam.save_fxn
CloudCovErr.decam.get_catnames
CloudCovErr.decam.proc_ccd
CloudCovErr.decam.proc_all
```
