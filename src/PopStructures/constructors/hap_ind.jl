# base types for individuals in age-structured populations
mutable struct HapIndividual <: AbstractIndividual
    age::Int
    nloci::Int                   # number of loci per haploid genome
    genome::Chromosome           # haploid genome
    genotype::Array{Float64,1}   # vector of genotypes
    phenotype::Array{Float64,1}  # vector of phenotypes
    fitness::Array{Float64,1}    # (survival, fertility)

    # constructors
    function HapIndividual(;age::Int=0, n::Int=1,
                        g::Vector{Float64}=Float64[],
                        p::Vector{Float64}=Float64[],
                        f::Vector{Float64}=[1.0, 1.0])
        i = new()
        i.age = age
        i.nloci = n
        i.genome = Chromosome(n)
        i.genotype = g
        i.phenotype = p
        i.fitness = f
        return i
    end
end
