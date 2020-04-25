# Simulation of population dynamics
#
# Jeremy Van Cleve <jeremy.vancleve@gmail.com>
# Daniel Priego Espinosa <dpriego87@gmail.com>

module JPop

using Distributions, LinearAlgebra

export
# interface
HaploidPop, DiploidPop, HaploDiploidPop,
# Basic types
Chromosome, HapIndividual, DipIndividual, Population,
# population quantities
mean_genotype, mean_phenotype, age_distribution,
# population dynamics
next_gen!, setInitFreq!, recombine

include("interface.jl")
include("utils.jl")
include("PopStructures/PopStructures.jl")
using .PopStructures
include("PopQuantities/PopQuantities.jl")
using .PopQuantities
include("EnvDynamics/EnvDynamics.jl")
using .EnvDynamics
include("PopDynamics/PopDynamics.jl")
using .PopDynamics

end # module
