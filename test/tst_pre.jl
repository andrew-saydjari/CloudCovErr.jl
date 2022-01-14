module tst_pre
    using Test
    using cloudCovErr

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

    end
end
