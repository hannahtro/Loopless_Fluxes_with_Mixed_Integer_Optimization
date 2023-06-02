# if success(`cmp --quiet model_loop_simple_model.lp model_vector_simple_model.lp`)
#     println("same")
# else
#     println("different")
# end

if success(`cmp --quiet model_loop.lp model_vector.lp`)
    println("same")
else
    println("different")
end

diff_output = read(ignorestatus(`diff model_loop.lp model_vector.lp`), String)
@show diff_output