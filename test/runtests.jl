using Test
using CloudCovErr

tests = [
    "tst_preprocess.jl",
    "tst_cov.jl",
    "tst_perstar.jl",
    "tst_decam.jl"
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
