<a name="readme-top"></a>

### README 

This README gives a short summary of the code used in my master's thesis "Mixed-Integer Optimization for Loopless Flux Distributions in Metabolic Networks".

### File Structure
The code for the different methods is located in `/src`. The code for FBA, ll-FBA and for the decomposition with combinatorial Benders' and with no-good cuts is on the main branch. 
To run ll-FBA on the model `e_coli_core`, type the following in the julia REPL:
```julia 
include("loopless_fba.jl")
loopless_fba_data("e_colicore")
```
To solve ll-FBA on the model `e_coli_core` with the combinatorial Benders' decomposition, run in the REPL:
```julia 
include("cuts_decomposition.jl")
combinatorial_benders_data("e_coli_core", fast=true, indicator=true, big_m=false)
```
If `fast=true`, we use combinatorial Benders' cuts, otherwise, no-good-cuts will be used. With `indicator=true` and `big_m=false`, the structure of the master problem is defined. If we want to block several cycles at once, we set `multiple_mis=0.5`, where we block `k` cycles and the number `k` is 0.5% of the number of reactions in the model.

We provide several tests in `/test` including tests for FBA, ll-FBA and for the decomposition with combinatorial Benders' and with no-good cuts. To run the tests, start a session in the folder with `julia --project` and run the command: `include("run_tests.jl")`.

The files to reproduce the experiments with the methods on the branch are in `/experiments`. The models used for the experiments are located in `/molecular_models`. The files for the thesis are found in `/latex`. Some code we did not use in the thesis is still in `/not_used`, including some experiments with MOMA.

### Branches
|  Branch | Content |
|---|---|
| enzyme_data | includes the randomized data generation for the enzyme models|
| disjunctive_programming | contains the code for the big-M reformulation and the convex-hull formulation using DisjunctiveProgramming.jl |
| intersection_cuts | contains the preliminary code for intersection cuts and the code for plots; the test file covers several tests including the implementation of an intersection cut for a simple MIP |
| plots_and_tables | contains the code for the performance plots and the code for selecting a subset of columns for the tables |
| cycle_free_flux | contains the code for blocking cycles and the code for generating the box plot | 
| constraint_handler | contains the preliminary experiments with using a constraint handler for the combinatorial Benders' approach|



