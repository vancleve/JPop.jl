module PopStructures
using JPop: AbstractIndividual,
            AbstractPopulation, HaploDiploidPop, HaploidPop, DiploidPop,
            rround
using Distributions
export Chromosome, HapIndividual, DipIndividual, Population

# include("./constructors/age_ind.jl")
include("./constructors/genome.jl")
include("./constructors/hap_ind.jl")
include("./constructors/dip_ind.jl")
include("./init_chr.jl")
include("./constructors/population.jl")

end
