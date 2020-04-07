# type for population, which is a collection of individuals,
# and population-level properties
mutable struct Population
    size::Int
    pheno_func!::Function            # map genotype to phenotype
    npheno::Int                   # number of phenotypes for each individual
    fit_func!::Function              # map phenotype to fitness
    mut_func!::Function              # mutate genotype
    members::Array{Individual,1}
    members_prev::Array{Individual,1}
    fitness::Array{Float64,2}
    mean_fit::Array{Float64,1}
    env_func!::Function
    env_state::Array{Float64,1}

    # Constructor
    ## geno_func function is used to initialize genotypes
    ## e.g., ()->[rand()] initializes each genotype to a random value in (0,1)
    function Population(size::Int,
                        pheno_func::Function, npheno::Int, fit_func::Function,
                        mut_func::Function, geno_func0::Function,
                        env_func::Function, env0::Array{Float64,1})

        members = Array{Individual}(undef,size)
        for i=1:size
            g = geno_func0(i)
            ind = Individual(0, g, Array{Float64}(undef,npheno))
            pheno_func(ind, env0)
            members[i] = ind
        end
        pop = new()
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
        return pop
    end
end
