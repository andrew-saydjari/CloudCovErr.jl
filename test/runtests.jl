using Test
using cloudCovErr

tests = [
    "tst_pre.jl",
    "tst_cov.jl",
    #"tst_plots.jl"
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
