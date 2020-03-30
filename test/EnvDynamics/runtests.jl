const envdir = dirname(@__FILE__)
tests = [
]

@testset "EnvDynamics" begin
    for t in tests
        tp = joinpath(envdir, "$(t).jl")
        include(tp)
    end
end
