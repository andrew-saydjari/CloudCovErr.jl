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
cloudCovErr.decam.read_decam
cloudCovErr.decam.read_crowdsource
cloudCovErr.decam.inject_rename
cloudCovErr.decam.load_psfmodel_cs
cloudCovErr.decam.save_fxn
cloudCovErr.decam.get_catnames
cloudCovErr.decam.proc_ccd
cloudCovErr.decam.proc_all
```
