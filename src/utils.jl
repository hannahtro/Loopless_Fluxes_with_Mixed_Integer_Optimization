# parse string to array, 
function parse_array_as_string(str)
    str = str[1]
    str = str[2:end-1]
    list = split(str, ", ")
    list = [parse(Float64, i) for i in list]
    # @show list
end 