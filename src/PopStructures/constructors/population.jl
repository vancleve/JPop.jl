# type for population, which is a collection of individuals,
# and population-level properties
mutable struct Population{T<:AbstractIndividual} <: AbstractPopulation
    size::Int
    pheno_func!::Function            # map genotype to phenotype
    npheno::Int                      # number of phenotypes for each individual
    fit_func!::Function              # map phenotype to fitness
    mut_func!::Function              # mutate genotype
    members::Array{T,1}
    members_prev::Array{T,1}
    fitness::Array{Float64,2}
    mean_fit::Array{Float64,1}
    env_func!::Function
    env_state::Array{Float64,1}
    rec_rf::Function        # function to set recombination rate vector
    rec_hat::Float64        # Average recombination rate
    mut_rf::Function        # function to set mutation rate vector
    mu_hat::Float64         # Average mutation rate
    # Constructor
    ## geno_func function is used to initialize genotypes
    ## e.g., ()->[rand()] initializes each genotype to a random value in (0,1)
    function Population(size::Int, ::Type{T},
                pheno_func::Function, npheno::Int,
                fit_func::Function, mut_func::Function,
                geno_func0::Function, env_func::Function,
                env0::Array{Float64,1}, nloci0::Int,
                rec_rf::Function, mut_rf::Function;
                # haploid fraction in a HaploDiploidPop
                hf::Float64=0.5,
                # average recombination rate
                rec_hat::Float64=0.5,
                # average mutation rate
                mu_hat::Float64=0.0) where {T<:AbstractPopulation}

        if T == HaploDiploidPop
            members = Array{AbstractIndividual}(undef,size)
            k = rround(size * hf)
            for i in 1:k
                g = geno_func0(i)
                ind = HapIndividual(age=0, n=nloci0, g=g,
                                p=Array{Float64}(undef,npheno))
                pheno_func(ind, env0)
                set_genome!(ind, rec_rf, rec_hat, mut_rf, mu_hat)
                members[i] = ind
            end
            for i in k+1:size
                g = geno_func0(i)
                ind = DipIndividual(age=0, n=nloci0, g=g,
                                p=Array{Float64}(undef,npheno))
                pheno_func(ind, env0)
                set_genome!(ind, rec_rf, rec_hat, mut_rf, mu_hat)
                members[i] = ind
            end
            S = AbstractIndividual
        else
            S = T == HaploidPop ? HapIndividual : DipIndividual
            members = Array{S}(undef,size)
            for i=1:size
                g = geno_func0(i)
                ind = S(age=0, n=nloci0, g=g,
                        p=Array{Float64}(undef,npheno))
                pheno_func(ind, env0)
                set_genome!(ind, rec_rf, rec_hat, mut_rf, mu_hat)
                members[i] = ind
            end
        end
        pop = new{S}()
        pop.size = size
        pop.pheno_func! = pheno_func
        pop.npheno = npheno
        pop.fit_func! = fit_func
        pop.mut_func! = mut_func
        pop.members = members
        pop.members_prev = copy(members)
        pop.fitness = hcat(fill(1, size), fill(1/size, size))
        pop.mean_fit = [1.0, 1/size]
        pop.env_func! = env_func
        pop.env_state = env0
        pop.rec_rf = rec_rf
        pop.rec_hat = rec_hat
        pop.mut_rf = mut_rf
        pop.mu_hat = mu_hat
        return pop
    end
end
