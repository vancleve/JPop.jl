# set initial frequency
function setInitFreq!(p, freq0, gvals = [1.0, 0.0])
    n = p.size
    k = rround(n * freq0)
    for i in 1:k
        p.members[i].genotype = [gvals[1]]
        for l in 1:p.members[i].nloci
            p.members[i].genome.loci_values[l] = gvals[1]
        end
    end
    for i in k+1:n
        p.members[i].genotype = [gvals[2]]
        for l in 1:p.members[i].nloci
            p.members[i].genome.loci_values[l] = gvals[2]
        end
    end
end
