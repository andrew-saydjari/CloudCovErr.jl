using Test
using cloudCovErr

tests = [
    "tst_preprocess.jl",
    "tst_cov.jl",
    "tst_perstar.jl"
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
