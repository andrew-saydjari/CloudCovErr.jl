# CloudCovErr.jl <img src="docs/src/assets/logo.png" alt="CloudCovErr Logo" width="100" align="right"/>

[![][docs-dev-img]][docs-dev-url]
[![][action-img]][action-url]
[![][codecov-img]][codecov-url]

Pipeline for debiasing and improving error bar estimates for photometry on top of structured/filamentary background. The procedure first estimates the covariance matrix of the residuals from a previous photometric model and then computes corrections to the estimated flux and flux uncertainties.

## Installation

**CloudCovErr** is a registered package so a stable version can be installed using `Pkg.add`.

```julia
import Pkg
Pkg.add("CloudCovErr")
```

For the most recent development version, install directly from the GitHub

```julia
import Pkg
Pkg.add("url=https://github.com/andrew-saydjari/CloudCovErr.jl")
```

Currently, we only support compatibility with linux and macOS in order to easily interface with dependencies of [crowdsource](https://github.com/schlafly/crowdsource). Due to older versions of Julia bundling outdated libstcd++, we only support Julia 1.6+ again to make interfacing with python-based photometric pipelines easier (see [issue](https://github.com/JuliaLang/julia/issues/34276)). However, workarounds exist for both problems. Please open an issue if there is some compatibility you would like supported.  

## Documentation

Detailed documentation can be found [here][docs-dev-url].

<!-- URLS -->
[action-img]: https://github.com/andrew-saydjari/CloudCovErr.jl/workflows/Unit%20test/badge.svg
[action-url]: https://github.com/andrew-saydjari/CloudCovErr.jl/actions

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://andrew-saydjari.github.io/CloudCovErr.jl/dev/

[codecov-img]: https://codecov.io/github/andrew-saydjari/CloudCovErr.jl/coverage.svg?branch=main
[codecov-url]: https://codecov.io/github/andrew-saydjari/CloudCovErr.jl?branch=main

## Contributing and Questions

This is a new piece of software. [Filing an
issue](https://github.com/andrew-saydjari/CloudCovErr.jl/issues/new) to report a
bug or request a feature is extremely valuable in helping us prioritize what to work on, so don't hesitate.
