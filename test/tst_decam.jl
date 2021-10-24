module decam_tests
    using Test
    using cloudCovErr
    using cloudCovErr.decam
    using Random

    @testset "imagePrep" begin
        ttt_testim = ones(33,33)
        ttt_maskim = zeros(Bool,33,33)
        ttt_maskim[16,16] = true
        ttt_bimage = ones(33,33)
        ttt_bimageI = ones(Int,33,33)
        ttt_testim2 = zeros(33,33)
        ttt_maskim2 = zeros(Bool,33,33)
        ttt_goodpix = ones(Bool,33,33)
        prelim_infill!(ttt_testim,ttt_maskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_maskim2, ttt_goodpix, widx = 19, widy=19)
        @test ttt_testim2[16,16] == 1.0

        ttt_testim2 = ones(2,2)
        ttt_maskim2 = ones(Bool,2,2)
        ttt_skyim3 = 10*ones(2,2)
        add_sky_noise!(ttt_testim2,ttt_maskim2,ttt_skyim3,4,seed=2021)
        @test ttt_testim2 == [0.25 2.25; 1.25 2.5]

        psfstamp = zeros(31,31)
        psfstamp[16,16] = 1
        ttt_maskim2 = zeros(Bool,51,51);
        gen_mask_staticPSF!(ttt_maskim2,psfstamp,[26],[26],[200])
        @test ttt_maskim2[26,26]

        # kpsf2d=zeros(Bool,3,3)
        # kpsf2d[2,2] = true
        # kpsf2d[2,3] = true
        # kstar = kpsf2d[:]
        # k = .!kstar
        # rng=MersenneTwister(2021)
        # data_in = randn(rng,3,3)
        # data_w = ones(3,3)
        # stars_in = 2*ones(3,3)
        # cov_loc = [
        # 47.2245   1.33946   0.628186   0.841306   0.288469  -0.437706   0.754434  -0.245601  -0.150857;
        #  0.0     47.2602    1.39354    0.570245   0.871201   0.243023  -0.457453   0.806491  -0.209147;
        #  0.0      0.0      47.166      1.40991    0.471935   0.767516   0.251835  -0.398556   0.772153;
        #  0.0      0.0       0.0       47.1832     1.44052    0.461676   0.757649   0.282371  -0.314066;
        #  0.0      0.0       0.0        0.0       47.2298     1.34296    0.471258   0.766785   0.218381;
        #  0.0      0.0       0.0        0.0        0.0       47.2845     1.36993    0.492573   0.715498;
        #  0.0      0.0       0.0        0.0        0.0        0.0       47.1459     1.3015     0.425096;
        #  0.0      0.0       0.0        0.0        0.0        0.0        0.0       47.1451     1.30432;
        #  0.0      0.0       0.0        0.0        0.0        0.0        0.0        0.0       47.2068;
        # ]
        # μ = [0.005037597380578518
        #      0.003629490500316024
        #      0.0026571564376354218
        #      0.0032599247060716152
        #      0.002772299339994788
        #      0.0034678715746849775
        #      0.0051109688356518745
        #      0.0037797277327626944
        #      0.003628455102443695]
        # psft = [
        # 0.0274769  0.0360258  0.032605;
        # 0.0330318  0.0403823  0.0339776;
        # 0.0301219  0.033388   0.0255433;
        # ]
        # @test all((condCovEst_wdiag(cov_loc,μ,kstar,kpsf2d,data_in,data_w,stars_in,psft) .- [133.725  132.769  5.59213e-5  -10.7856  0.906103  4.94119e-5]) .< 1e-2)

        ttt_residimIn = ones(51,51)
        ttt_maskim = zeros(Bool,51,51)
        ttt_maskim[26,26] = true
        ttt_w_im = ones(51,51)
        ttt_mod_im = ones(51,51)
        ttt_skyim = ones(51,51);
        data_in, data_w, stars_in, kmasked2d = stamp_cutter(26.1,26.1,ttt_residimIn,ttt_w_im,ttt_mod_im,ttt_skyim,ttt_maskim)
        @test kmasked2d[17,17]
    end
end
