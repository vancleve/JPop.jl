struct Chromosome
    loci_values::Vector{Float64} # vector of allele values
    loci_ids::Vector{String}     # vector of identifiers
    rec_rate::Vector{Float64}    # vector of recombination rates
    mu::Vector{Float64}     # vector of mutation rates, these three vectors have the same length
    coop_loci::BitArray # vector of boolean flags to mark cooperation related loci
    Chromosome(n::Int) = new(Vector{Float64}(undef,n),
                            Vector{String}(undef,n),
                            Vector{Float64}(undef,n),
                            Vector{Float64}(undef,n),
                            falses(n))
    function Chromosome(loci_values::Vector{Float64},
                        loci_ids::Vector{String},
                        rec_rate::Vector{Float64},
                        mu::Vector{Float64},
                        coop_loci::AbstractArray{Bool})
        if (size(loci_values) == size(loci_ids) == size(rec_rate) == size(mu) == size(coop_loci))
            new(loci_values, loci_ids, rec_rate, mu, coop_loci)
        else
            throw("All fields must have the same size")
        end
    end
end
