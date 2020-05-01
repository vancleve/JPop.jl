function set_genome!(i::HapIndividual, rec_rf::Function, rec_hat::Float64, mut_rf::Function, mu_hat::Float64; tag="A")
    n = i.nloci
    chr = i.genome
    for ni=1:n
        chr.loci_ids[ni] = "$tag$ni"
    end
    chr.rec_rate .= rec_rf(n,rec_hat)
    chr.mu .= mut_rf(n,mu_hat)
    nothing
end

function set_genome!(i::DipIndividual, rec_rf::Function, rec_hat::Float64, mut_rf::Function, mu_hat::Float64; tag="A")
    n = i.nloci
    for nhap=1:2
        chr = i.genome[nhap]
        for ni=1:n
            i.chr.loci_ids[ni] = "$tag$ni"
        end
        chr.rec_rate .= rec_rf(n,rec_hat)
        chr.mu .= mut_rf(n,mu_hat)
    end
    nothing
end
