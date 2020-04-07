@testset "Population" begin
    # population size
    size = 1000
    # number of phenotypes
    npheno = 1
    # genotype->phenotype map,
    function phenof!(i::Individual,e)
        i.phenotype[1] = 2.0
    end
    # fitness function for Wright-Fisher model
    function fitf!(i::Individual,e)
        # survival: non-overlapping generations (at each iteration, previous
        # generation individuals die and are replaced by their descendents)
        i.fitness[1] = 0.0
        # fertility set constant
        i.fitness[2] = 1.0
        nothing
    end
    # mutation function
    function mutf!(offspring::Individual, parent::Individual)
        copy!(offspring.genotype, parent.genotype)
        nothing
    end
    # initial genotype function
    g0 = (i) -> [0.0]
    # env update function
    envf = (e) -> [1.0]
    # initial env state
    env = [1.0]
    pop = Population(size, phenof!, npheno, fitf!, mutf!, g0, envf, env0);
    ni = rand(1:size)
    @test pop.size == size
    m = pop.members[ni]
    @test m.phenotype[ni] == phenof!(m,0)
    @test m.genotype == g0(m)
end
