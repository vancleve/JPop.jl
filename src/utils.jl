# round to nearest integer randomly
function rround(x)
    return floor(Int,x) + (rand() < x - floor(Int,x))
end
