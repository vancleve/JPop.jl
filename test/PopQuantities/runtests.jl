const quantdir = dirname(@__FILE__)
tests = [
]

@testset "PopQuantities" begin
    for t in tests
        tp = joinpath(quantdir, "$(t).jl")
        include(tp)
    end
end
