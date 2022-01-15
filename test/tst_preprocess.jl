module tst_pre
    using Test
    using CloudCovErr
    using Random

    # Random generators do not promise the same streams between
    # different Julia releases... that is highly frustrating.

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
        @test out[2] == 12:17
        @test length(out[3]) == 20
        out = im_subrng(1,1,cx,cy,16,16,2,2,4,4,1,1,4,4)
        @test out[1] == 0:5
        out = im_subrng(4,4,cx,cy,15,15,2,2,4,4,1,1,4,4)
        @test out[1] == 12:16

        refin = [
         11.0  12.0   11.0 ;
         32.0  60.0   19.0 ;
         22.0  25.0   20.0 ;
         ]
        add_noise!(refin,2;seed=2021)
        ref = [
         12.0   9.5  10.0 ;
         31.5  70.5  18.0 ;
         19.0  25.5  19.5 ;
         ]
        @test abs(refin[2,2]-70.5) < 10

        skyim = 10 .*ones(3,3)
        maskim = zeros(Bool,3,3)
        maskim[1,:] .= true
        testim2 = [
          1.0  1.0  1.0;
         -1.0  1.0  1.0;
          1.0  1.0  1.0;
        ]
        add_sky_noise!(testim2,maskim,skyim,4;seed=2021)
        ref = [
         0.25  1.25  2.25 ;
         1.0  1.0   -1.0 ;
         -1.0  1.0   -1.0 ;
         ]
        @test all(abs.(ref) .< 5)

        ttt_testim = ones(33,33)
        ttt_bmaskim = zeros(Bool,33,33)
        ttt_bmaskim[16,16] = true
        ttt_bimage = ones(33,33)
        ttt_bimageI = ones(Int,33,33)
        ttt_testim2 = zeros(33,33)
        ttt_bmaskim2 = zeros(Bool,33,33)
        ttt_goodpix = ones(Bool,33,33)
        prelim_infill!(ttt_testim,ttt_bmaskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_bmaskim2,ttt_goodpix,"N4",widx=19,widy=19)
        @test ttt_testim2[16,16] == 1.0
        prelim_infill!(ttt_testim,ttt_bmaskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_bmaskim2,ttt_goodpix,"N4",widx=1,widy=1)
        @test ttt_testim2[16,16] == 1.0
        prelim_infill!(ttt_testim,ttt_bmaskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_bmaskim2,ttt_goodpix,"N4",widx=19,widy=19,ftype=64)
        @test ttt_testim2[16,16] == 1.0

        psfstamp = zeros(31,31)
        psfstamp[16,16] = 1
        ttt_maskim2 = zeros(Bool,51,51);
        gen_mask_staticPSF!(ttt_maskim2,psfstamp,[26],[26],[200])
        @test ttt_maskim2[26,26]

        psfstamp = zeros(33,33)
        psfstamp[16,16] = 100
        psfstamp[16,1] = 2

        psfstamp1 = zeros(31,31)
        psfstamp1[16,16] = 100
        psfstamp1[1,8] = 1

        ttt_maskim2 = zeros(Bool,51,51);
        gen_mask_staticPSF2!(ttt_maskim2,psfstamp,psfstamp1,[26],[26],[20])
        @test ttt_maskim2[26,26]
        @test !ttt_maskim2[26,26-16]

        psfstamp = zeros(33,33)
        psfstamp[17,17] = 100
        psfstamp[17,1] = 2

        psfstamp1 = zeros(31,31)
        psfstamp1[16,16] = 100
        psfstamp1[1,8] = 1

        ttt_maskim2 = zeros(Bool,51,51);
        gen_mask_staticPSF2!(ttt_maskim2,psfstamp,psfstamp1,[26],[26],[200])
        @test ttt_maskim2[26,26]
        @test ttt_maskim2[26,26-16]
    end
end
