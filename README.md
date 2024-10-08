## Appraisal + iterative optimisation of mosquito surveillance in Australia

Welcome to the repo for Chapter 6 of my thesis, *Appraisal and iterative selection 
of sites for surveillance of Japanese encephalitis virus in mosquitoes*. 

In `/code/` directory you will find:

- `main.R`: opens up libraries and reads in data

- `genetic_algot.R`: includes all functions for genetic algorithm, including 
some plotting and analysis functions

- `*_setup.R`: organises data for specific jurisdictions, one of Victoria (`vic`) and 
Western Australia (`wa`)

- `*_pool.R`, `*_neigh.R`, `*_pareto.R`, `*_greedy.R`, etc.: runs the genetic algorithm 
for a particular jurisdiction, for a particular set of settings (e.g., `pool` is for 
scripts with different pool sizes); saves results (in `/output/`) to `.csv`s of area 
under the Pareto front and `.rds` objects of simulation times and members of the final 
simulated Pareto fronts

- `*_outs.R`: generates figures for the body of the chapter for a particular jurisdiction. 
(`vic_outs.R` also includes the code for a couple of supp figures investigating the area 
under the Pareto front metric)

- `sa4_appraisal.R`: completes the SA4-level "appraisal" of current mosquito surveillance effort
and creates the map and table outputs for the chapter

(Don't be afraid of the `trapezoid`: these snippets/results are more up to date as 
I changed the method I was using to calculate the area under the front)
