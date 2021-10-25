module tst_plots
    using Test
    using cloudCovErr


    @testset "Plots" begin
        cov_out = ones(2,2)
        @test cov_as_img(cov_out) != 1
        @test plot_cov_compare(cov_out, cov_out) != 1
    end
end
