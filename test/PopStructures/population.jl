@testset "Population" begin
    # population size
    size = 1000
    # number of phenotypes
    npheno = 1
    # genotype->phenotype map,
    function phenof!(i::HapIndividual,e)
        i.phenotype[1] = 2.0
    end
    # fitness function for Wright-Fisher model
    function fitf!(i::HapIndividual,e)
        # survival: non-overlapping generations (at each iteration, previous
        # generation individuals die and are replaced by their descendents)
        i.fitness[1] = 0.0
        # fertility set constant
        i.fitness[2] = 1.0
        nothing
    end
    # mutation function
    function mutf!(offspring::HapIndividual, parent::HapIndividual)
        copy!(offspring.genotype, parent.genotype)
        nothing
    end
    # initial genotype function
    g0 = (i) -> [0.0]
    # env update function
    envf = (e) -> [1.0]
    # initial env state
    env0 = [1.0]
    nloci = 1

    function recf(n::Int, rec_hat::Float64)
        rand(Normal(rec_hat,0.25*rec_hat),n)
    end

    function mut_rf(n::Int, mu_hat::Float64)
        rand(Normal(mu_hat,0.25*mu_hat),n)
    end

    pop = Population(size, HaploidPop, phenof!, npheno, fitf!, mutf!, g0,
                    envf, env0, nloci, recf, mut_rf);
    ni = rand(1:size)
    @test pop.size == size
    m = pop.members[ni]
    @test m.phenotype[1] == phenof!(m,0)
    @test m.genotype == g0(m)
end
