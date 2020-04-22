# set initial frequency
function setInitFreq(p, freq0, gvals = [1.0, 0.0])
    n = p.size
    k = rround(n * freq0)
    for i in 1:k
        p.members[i].genotype = [gvals[1]]
    end
    for i in k+1:n
        p.members[i].genotype = [gvals[2]]
    end
end
