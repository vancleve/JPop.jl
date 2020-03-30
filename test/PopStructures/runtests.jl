const structdir = dirname(@__FILE__)
tests = ["age_ind"
]

@testset "PopStructures" begin
    for t in tests
        tp = joinpath(structdir, "$(t).jl")
        include(tp)
    end
end
