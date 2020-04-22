
###
### Main life cycle function:
### iterate through the lifecycle of the population once
### 1. fitness calculated
### 2. survival
### 3. fertility
### 4. update environmental state
###
function next_gen!(pop::Population)
    # save initial population state in "prev" vector
    copy!(pop.members_prev, pop.members)

    # update fitness values
    calc_pheno_fitness!(pop)

    # get normalized fertility
    norm_fert = pop.fitness[:,2] / (pop.mean_fit[2] * pop.size)

    # set categorical distribution using normalized fertilities ("sampler" uses "AliasTable")
    fertdist = sampler(Categorical(norm_fert))

    # survival and reproduction
    for i=1:pop.size
        if pop.members[i].fitness[1] > rand()
            # individual survives and ages
            pop.members[i].age += 1
        else
            # individual dies and is replaced by random new born
            pop.members[i].age = 0
            parent = rand(fertdist)
            pop.mut_func!(pop.members[i], pop.members_prev[parent]) # offspring, parent
        end
    end

    # update environmental state
    update_env_state!(pop)

end
