# base types for individuals in age-structured populations
mutable struct Individual
    age::Int
    genotype::Array{Float64,1}   # vector of genotypes
    phenotype::Array{Float64,1}  # vector of phenotypes
    fitness::Array{Float64,1}    # (survival, fertility)

    # constructors
    Individual() = new(0, Float64[], Float64[], [1.0, 1.0])
    Individual(age::Int) = new(age, Float64[], Float64[], [1.0, 1.0])
    Individual(age::Int, genotype::Vector{Float64}) = new(age, genotype, Float64[], [1.0, 1.0])
    Individual(age::Int, genotype::Vector{Float64}, phenotype::Vector{Float64}) = new(age, genotype, phenotype, [1.0, 1.0])
    Individual(age::Int, genotype::Vector{Float64}, phenotype::Vector{Float64}, fitness::Vector{Float64}) = new(age, genotype, phenotype, fitness)
end
