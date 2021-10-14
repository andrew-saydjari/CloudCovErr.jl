using Test
using disCovErr

tests = [
    "tst_utils.jl",
]

for t in tests
    @testset "$t" begin
        include(t)
    end
end
