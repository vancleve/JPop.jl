using JPop
using Test, Distributions

const testdir = dirname(@__FILE__)

tests = [
    "EnvDynamics",
    "PopDynamics",
    "PopQuantities",
    "PopStructures"
]

@testset "JPop.jl" begin
    for t in tests
        tp = joinpath(testdir, "$(t)/runtests.jl")
        include(tp)
    end
end
