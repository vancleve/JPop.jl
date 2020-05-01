module PopDynamics
# import base functions for multiple dispatch
import Base.copy, Base.copy!
using JPop: rround, AbstractIndividual
using JPop.PopStructures: Chromosome, HapIndividual, DipIndividual, Population
using JPop.EnvDynamics: update_env_state!
using JPop.PopQuantities: calc_pheno_fitness!

# using ..EnvDynamics: update_env_state!
# using ..PopQuantities: calc_pheno_fitness!
# using ..PopStructures: Individual, Population
using Distributions
export next_gen!, setInitFreq!, recombine

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

function recombine(i::HapIndividual, j::HapIndividual)
    parent_chr = Vector{Chromosome}(undef,2)
    parent_chr[1] = i.genome
    parent_chr[2] = j.genome
    r_ind = 1
    r_count = 0             # Number of recombination events per reproduction
    chosen_chr = parent_chr[r_ind]
    rec_g = Chromosome(i.nloci)
    for l=1:i.nloci
        if rand() < parent_chr[r_ind].rec_rate[l]
            r_count += 1
            r_ind = 1 + r_count % 2
            chosen_chr = parent_chr[r_ind]
        end
        rec_g.loci_values[l] = chosen_chr.loci_values[l]
        rec_g.loci_ids[l] = chosen_chr.loci_ids[l]
        rec_g.rec_rate[l] = chosen_chr.rec_rate[l]
        rec_g.mu[l] = chosen_chr.mu[l]
        rec_g.coop_loci[l] = chosen_chr.coop_loci[l]
    end
    return r_count, rec_g
end

function recombine(i::DipIndividual, j::DipIndividual)
    parents = Vector{DipIndividual}(undef,2)
    parents[1] = i
    parents[2] = j
    parent_chr = Vector{Chromosome}(undef,2)
    r_count = 0             # Number of recombination events per reproduction
    rec_dip_genome = Vector{Chromosome}(undef,2)
    for pn in 1:2
        p = parents[pn]
        parent_chr[1] = p.genome[1]
        parent_chr[2] = p.genome[2]
        r_ind = 1
        chosen_chr = parent_chr[r_ind]
        rec_hap_g = Chromosome(p.nloci)
        for l=1:p.nloci
            if rand() < parent_chr[r_ind].rec_rate[l]
                r_count += 1
                r_ind = 1 + r_count % 2
                chosen_chr = parent_chr[r_ind]
            end
            rec_hap_g.loci_values[l] = chosen_chr.loci_values[l]
            rec_hap_g.loci_ids[l] = chosen_chr.loci_ids[l]
            rec_hap_g.rec_rate[l] = chosen_chr.rec_rate[l]
            rec_hap_g.mu[l] = chosen_chr.mu[l]
            rec_hap_g.coop_loci[l] = chosen_chr.coop_loci[l]
        end
        rec_dip_genome[pn] = rec_hap_g
    end
    return r_count, rec_dip_genome
end

include("./initializers.jl")
include("./updaters.jl")

end
