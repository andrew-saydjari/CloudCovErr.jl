module tst_pre
    using Test
    using cloudCovErr
    using Random

    @testset "imagePrep" begin
        out = kstar_circle_mask(3;rlim=1)
        ref = [
         true   false   true   ;
         false  false   false  ;
         true   false   true   ;
         ]
        @test out == ref

        ref = [
         1.1   1.2   1.1 ;
         3.2  6.0   1.9  ;
         2.2   2.5   2.0 ;
         ]
        @test findmaxpsf(ref;thr=20) == 6.25

        cx = (0:17) * ones(18)'
        cy = ones(18) * (0:17)'
        out = im_subrng(2,4,cx,cy,16,16,2,2,4,4,1,1,4,4)
        @test out[1] == 4:9
        @test out[1] == 12:17
        @test length(out[3]) == 20

        refin = [
         11.0  12.0   11.0 ;
         32.0  60.0   19.0 ;
         22.0  25.0   20.0 ;
         ]
        add_noise!(refin,2;seed=2021)
        ref = [
         10.5   9.5  14.0 ;
         26.0  64.5  18.0 ;
         21.0  21.0  17.5 ;
         ]
        @test ref == refin

        rng = MersenneTwister(2022)
        skyim = 10 .*ones(3,3)
        maskim = zeros(Bool,3,3)
        maskim[1,:] .= true
        testim2 = rand(rng,[-1.0,1.0],3,3)
        add_sky_noise!(testim2,maskim,skyim,4;seed=2021)
        ref = [
         1.5  3.25  1.5 ;
         -1.0  1.0   1.0 ;
         1.0  1.0   1.0 ;
         ]
        @test testim2 == ref

    end
end
