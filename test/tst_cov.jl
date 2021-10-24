module tst_cov
    using Test
    using cloudCovErr

    @testset "CovEst" begin
        X, Y = cov_construct(ones(256,256),[128],[128])
        @test (X == zeros(1,33^2,33^2)).&(Y==ones(1,33^2))

        Np = 5
        Δx = 2
        Δy = 2
        img = collect(reshape(1:400,20,20))
        (sx, sy) = size(img)
        @test (cov_construct(img, [1], [10])[1][1,1,1] - 11982.952403846151) < 1e-7
        @test (cov_construct(img, [20], [10])[1][1,1,1] - 11982.952403846151) < 1e-7
        @test (cov_construct(img, [25], [10])[1][1,1,1] - 11981.526682692307) < 1e-7
        @test (cov_construct(img, [-5], [10])[1][1,1,1] - 11983.214182692302) < 1e-7
    end
end
