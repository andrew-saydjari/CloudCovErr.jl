using Test
using cloudCovErr

tests = [
    "tst_cov.jl",
    "tst_decam.jl"
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
