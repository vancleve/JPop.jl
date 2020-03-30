const dyndir = dirname(@__FILE__)
tests = [
]

@testset "PopDynamics" begin
    for t in tests
        tp = joinpath(dyndir, "$(t).jl")
        include(tp)
    end
end
