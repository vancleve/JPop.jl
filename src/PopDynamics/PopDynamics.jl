module PopDynamics
# import base functions for multiple dispatch
import Base.copy, Base.copy!
using JPop: rround, update_env_state!, calc_pheno_fitness!, Chromosome, AbstractIndividual, Population
# using ..EnvDynamics: update_env_state!
# using ..PopQuantities: calc_pheno_fitness!
# using ..PopStructures: Individual, Population
using Distributions
export next_gen!, setInitFreq

# copy methods for chromosomes
function copy(i::Chromosome)
    return Chromosome(i.loci_values,i.loci_ids,i.rec_rate,i.mu,i.coop_loci)
end
# copy method for array of Chromosomes
function copy(x::Array{Chromosome,1})
    y = Array{Chromosome}(undef,size(x))
    for i=1:length(x)
        y[i] = copy(x[i])
    end
    return y
end

# copy method for Individual
function copy(i::T) where  {T<:AbstractIndividual}
    j = T(age=i.age, n=i.nloci, g=copy(i.genotype),
        p=copy(i.phenotype),
        f=copy(i.fitness))
    j.genome = copy(i.genome)
    return j
end
# copy inplace method for Individual
function copy!(i::T, j::T) where {T<:AbstractIndividual}
    i.age = j.age
    i.genome = copy(i.genome)
    copy!(i.genotype, j.genotype)
    copy!(i.phenotype, j.phenotype)
    copy!(i.fitness, j.fitness)
end

# copy method for Array of Individuals
function copy(x::Array{T,1}) where {T<:AbstractIndividual}
    y = Array{T}(undef,size(x))
    for i=1:length(x)
        y[i] = copy(x[i])
    end
    return y
end

# copy inplace method for Array of Individuals
function copy!(x::Array{T,1}, y::Array{T,1}) where {T<:AbstractIndividual}
    if (length(x) != length(y))
        error("Arrays must have equal length")
    end

    for i=1:length(x)
        copy!(x[i], y[i])
    end
end

include("./initializers.jl")
include("./updaters.jl")

end
