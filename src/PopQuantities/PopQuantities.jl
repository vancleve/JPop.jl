module PopQuantities
using ..PopStructures: Population
using Distributions
export mean_genotype, mean_phenotype, age_distribution, calc_pheno_fitness!

function calc_pheno_fitness!(pop::Population)
    pop.mean_fit = [0.0, 0.0]

    @inbounds for i = 1:pop.size
        pop.pheno_func!(pop.members[i], pop.env_state)
        pop.fit_func!(pop.members[i], pop.env_state)
        pop.fitness[i,1] = pop.members[i].fitness[1]
        pop.fitness[i,2] = pop.members[i].fitness[2]
        pop.mean_fit[1] += pop.members[i].fitness[1] / pop.size
        pop.mean_fit[2] += pop.members[i].fitness[2] / pop.size
    end
    nothing
end

function mean_genotype(pop::Population)
    ngeno = length(pop.members[1].genotype)
    pop_geno = Array{Float64}(undef,ngeno, pop.size)

    for i=1:pop.size
        pop_geno[:,i] = pop.members[i].genotype
    end

    return mean(pop_geno, dims=2)
end

function mean_phenotype(pop::Population)
    npheno = length(pop.members[1].phenotype)
    pop_pheno = Array{Float64}(undef,npheno, pop.size)

    for i=1:pop.size
        pop_pheno[:,i] = pop.members[i].phenotype
    end

    return mean(pop_pheno, 2)
end

function age_distribution(pop::Population)
    ages = Array{Float64}(undef,pop.size)

    for i=1:pop.size
        ages[i] = pop.members[i].age
    end
    return hist(ages, -0.5:1:maximum(ages)+0.5)
end

end
