module tst_cov
    using Test
    using cloudCovErr

    @testset "DECam" begin
        ref = "/src/data/decapsi/c4d_170119_085651_ood_r_v1.I.fits.fz"
        @test ref == cloudCovErr.inject_rename("/src/data/decaps/c4d_170119_085651_ood_r_v1.fits.fz")

        
    end
end
