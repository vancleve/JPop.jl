@testset "HapIndividual" begin
    i = HapIndividual()
    @test i.age == 0
    @test i.nloci == 1
    @test i.genotype == Float64[]
    @test i.phenotype == Float64[]
    @test i.fitness == [1.0,1.0]
    age = 17
    i = HapIndividual(age=age)
    @test i.age == age
    @test i.genotype == Float64[]
    @test i.phenotype == Float64[]
    @test i.fitness == [1.0,1.0]
    genotype = [0.25, 0.5]
    i = HapIndividual(age=age, g=genotype)
    @test i.age == age
    @test i.genotype == genotype
    @test i.phenotype == Float64[]
    @test i.fitness == [1.0,1.0]
    phenotype = reverse(genotype)
    i = HapIndividual(age=age, g=genotype, p=phenotype)
    @test i.age == age
    @test i.genotype == genotype
    @test i.phenotype == phenotype
    @test i.fitness == [1.0,1.0]
    fitness = [2.0 0; 0.0 0.5]*phenotype
    i = HapIndividual(age=age, g=genotype, p=phenotype, f=fitness)
    @test i.age == age
    @test i.genotype == genotype
    @test i.phenotype == phenotype
    @test i.fitness == fitness
end
