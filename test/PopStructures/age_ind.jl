@testset "Individual" begin
    i = Individual()
    @test i.age == 0
    @test i.genotype == Float64[]
    @test i.phenotype == Float64[]
    @test i.fitness == [1.0,1.0]
    age = 17
    i = Individual(age)
    @test i.age == age
    @test i.genotype == Float64[]
    @test i.phenotype == Float64[]
    @test i.fitness == [1.0,1.0]
    genotype = [0.25, 0.5]
    i = Individual(age, genotype)
    @test i.age == age
    @test i.genotype == genotype
    @test i.phenotype == Float64[]
    @test i.fitness == [1.0,1.0]
    phenotype = reverse(genotype)
    i = Individual(age, genotype,phenotype)
    @test i.age == age
    @test i.genotype == genotype
    @test i.phenotype == phenotype
    @test i.fitness == [1.0,1.0]
    fitness = [2.0 0; 0.0 0.5]*phenotype
    i = Individual(age, genotype,phenotype,fitness)
    @test i.age == age
    @test i.genotype == genotype
    @test i.phenotype == phenotype
    @test i.fitness == fitness
end
