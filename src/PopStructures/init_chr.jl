function set_genome!(i::HapIndividual, rec_rf::Function, rec_hat::Float64, mut_rf::Function, mu_hat::Float64; tag="A")
    n = i.nloci
    chr = Chromosome(n)
    for ni=1:n
        chr.loci_ids[ni] = "$tag$ni"
    end
    chr.rec_rate .= rec_rf(n,rec_hat)
    chr.mu .= mut_rf(n,mu_hat)
    i.genome = chr
    nothing
end

function set_genome!(i::DipIndividual, rec_rf::Function, rec_hat::Float64, mut_rf::Function, mu_hat::Float64; tag="A")
    n = i.nloci
    for nhap=1:2
        chr = Chromosome(n)
        for ni=1:n
            chr.loci_ids[ni] = "$tag$ni"
        end
        chr.rec_rate .= rec_rf(n,rec_hat)
        chr.mu .= mut_rf(n,mu_hat)
        i.genome[nhap] = chr
    end
    nothing
end
