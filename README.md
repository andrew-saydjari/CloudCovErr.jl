# cloudCovErr.jl <img src="docs/src/assets/logo.png" alt="cloudCovErr Logo" width="100" align="right"/>

[![][action-img]][action-url]
[![][docs-dev-img]][docs-dev-url]
[![][codecov-img]][codecov-url]

Pipeline for debiasing and improving error bar estimates for photometry on top of structured/filamentary background. The procedure first estimates the covariance matrix of the residuals from a previous photometric model and then computes corrections to the estimated flux and flux uncertainties.

## Installation

**cloudCovErr** is a registered package so it can be installed using `Pkg.add`.

```julia
import Pkg
Pkg.add("cloudCovErr")
```

## Documentation

Detailed documentation can be found [here][docs-dev-url].

<!-- URLS -->
[action-img]: https://github.com/andrew-saydjari/cloudCovErr.jl/workflows/Unit%20test/badge.svg
[action-url]: https://github.com/andrew-saydjari/cloudCovErr.jl/actions

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://andrew-saydjari.github.io/cloudCovErr.jl/dev/

[codecov-img]: https://codecov.io/github/andrew-saydjari/cloudCovErr.jl/coverage.svg?branch=main
[codecov-url]: https://codecov.io/github/andrew-saydjari/cloudCovErr.jl?branch=main

## Contributing and Questions

This is a new piece of software. [Filing an
issue](https://github.com/andrew-saydjari/cloudCovErr.jl/issues/new) to report a
bug or request a feature is extremely valuable in helping us prioritize what to work on, so don't hesitate.
