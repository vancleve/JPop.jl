mutable struct DipIndividual <: AbstractIndividual
    age::Int
    nloci::Int                   # number of loci per haploid genome
    genome::Vector{Chromosome}   # Vector of haploid genomes
    genotype::Array{Float64,1}   # vector of genotypes
    phenotype::Array{Float64,1}  # vector of phenotypes
    fitness::Array{Float64,1}    # (survival, fertility)

    # constructors
    function DipIndividual(;age::Int=0, n::Int=1,
                        g::Vector{Float64}=Float64[],
                        p::Vector{Float64}=Float64[],
                        f::Vector{Float64}=[1.0, 1.0])
        i = new()
        i.age = age
        i.nloci = n
        i.genome = Vector{Chromosome}(undef,2)
        i.genotype = g
        i.phenotype = p
        i.fitness = f
        return i
    end
end
