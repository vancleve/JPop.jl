# Simulation of population dynamics
#
# Jeremy Van Cleve <jeremy.vancleve@gmail.com>
# Daniel Priego Espinosa <dpriego87@gmail.com>

module JPop

using Distributions, LinearAlgebra

export
# Basic types
Individual, Population,
# population quantities
mean_genotype, mean_phenotype, age_distribution,
# population dynamics
next_gen!

include("PopStructures/PopStructures.jl")
using .PopStructures
include("PopQuantities/PopQuantities.jl")
using .PopQuantities
include("EnvDynamics/EnvDynamics.jl")
using .EnvDynamics
include("PopDynamics/PopDynamics.jl")
using .PopDynamics

end # module
