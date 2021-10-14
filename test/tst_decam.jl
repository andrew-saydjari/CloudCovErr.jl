module decam_tests
    using Test
    using disCovErr

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
        @test ttt_testim2 == [1.5 1.5; 3.25 2.5]

    end
end
