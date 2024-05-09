
function presolve_write_from_file(problem_input::String, problem_postsolve::String, reduced_problem::String)
    @assert isfile(problem_input)
    SCIP_PaPILO_jll.papilo() do exe
        run(`$exe presolve -f $problem_input -v $problem_postsolve -r $reduced_problem`)
    end
end