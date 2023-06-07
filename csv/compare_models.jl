# if success(`cmp --quiet model_loop_simple_model.lp model_vector_simple_model.lp`)
#     println("same")
# else
#     println("different")
# end

if success(`cmp --quiet model_vector_iAF692_laptop.lp model_vector_iAF692_server.lp`)
    println("same")
else
    println("different")
end

diff_output = read(ignorestatus(`diff model_vector_iAF692_laptop.lp model_vector_iAF692_server.lp`), String)
@show diff_output