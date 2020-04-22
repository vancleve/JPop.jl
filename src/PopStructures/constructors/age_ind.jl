# base types for individuals in age-structured populations
mutable struct HapIndividual <: AbstractIndividual
    age::Int
    genome::Chromosome
    genotype::Array{Float64,1}   # vector of genotypes
    phenotype::Array{Float64,1}  # vector of phenotypes
    fitness::Array{Float64,1}    # (survival, fertility)

    # constructors
    function HapIndividual(age::Int,
                        # function to initialize chromosome
                        chr_init::Function,
                        # Number of loci in a chromosome
                        n::Int,
                        genotype::Vector{Float64},
                        phenotype::Vector{Float64},
                        fitness::Vector{Float64})
        i = new()
        i.age = age
        i.genome = chr_init(n)
        i.genotype = genotype
        i.phenotype = phenotype
        i.fitness = fitness
        return i
    end
end

HapIndividual() = HapIndividual(0, chr_init, Float64[], Float64[], [1.0, 1.0])
HapIndividual(age) = HapIndividual(age, chr_init, Float64[], Float64[], [1.0, 1.0])
HapIndividual(age,chr_init, genotype) = HapIndividual(age, chr_init, genotype, Float64[], [1.0, 1.0])
HapIndividual(age,chr_init, genotype, phenotype) = HapIndividual(age, chr_init, genotype, phenotype, [1.0, 1.0])
HapIndividual(age,chr_init, genotype, phenotype, fitness) = HapIndividual(age, chr_init, genotype, phenotype, fitness)
