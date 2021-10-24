module utils_tests
    using Test
    using cloudCovErr

    @testset "CovEst" begin
        X, Y = cov_construct(ones(256,256),[128],[128])
        @test (X == zeros(1,33^2,33^2)).&(Y==ones(1,33^2))

    end
end
