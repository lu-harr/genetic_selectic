plot(vic_objective$potent, vic_objective$hpop)

ncell(vic_objective)
choose(ncell(vic_objective), 50) # how many sites in original appraisal for 2023-24? Or 2022-23?
choose(ncell(toy_objective), 5)

# I guess compare performances by plotting fronts next to each other?
# Stop when there hasn't been a change for 10 iters?