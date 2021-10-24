using Test
using cloudCovErr

tests = [
    "tst_utils.jl",
    "tst_decam.jl"
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
