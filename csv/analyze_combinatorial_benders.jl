using DataFrames
using CSV
using PyPlot

function plot_dual_bound(organism, file_name; csv=true)
    df = first(CSV.read(organism * "_" * file_name * ".csv", DataFrame),1)
    @show df

    objective_value = df[!,:objective_value]
    objective_values = df[!,:objective_values][1]
    objective_values = replace(objective_values, "Any[" => "")
    objective_values = replace(objective_values, "]" => "")
    objective_values = replace(objective_values, "," => "")
    objective_values = parse.(Float64, split(objective_values))
    dual_bounds = df[!,:dual_bounds][1]
    dual_bounds = replace(dual_bounds, "Any[" => "")
    dual_bounds = replace(dual_bounds, "]" => "")
    dual_bounds = replace(dual_bounds, "," => "")
    dual_bounds = parse.(Float64, split(dual_bounds))
    solution = df[!,:solution]
    time = df[!,:time]
    termination = df[!,:termination]
    iter = df[!,:iter][1] + 1

    # @show objective_value
    @show isapprox(objective_values[end],0)
    @show length(objective_values)
    @show objective_values[end]
    @show iter

    fig = plt.figure(figsize=(6.5, 3.5))
    ax = fig.add_subplot(111)
    # ax.plot(collect(1:length(dual_bounds)), dual_bounds)
    ax.plot(collect(1:length(objective_values)), objective_values)

    ylabel("dual bound")
    #locator_params(axis="y", nbins=4)
    xlabel("iteration")
    ax.grid()
    fig.tight_layout()
    f = matplotlib.ticker.FormatStrFormatter("%1.1f") # Define format of tick labels
    ax.xaxis.set_major_formatter(f) # Set format of tick labels

    if csv
        file = "../plots/" * organism * "_" * file_name * ".pdf"
        savefig(file)
    end
end 

organism = "iAF692"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name, csv=false)

# organism = "iAF692"
# file_name = "combinatorial_benders_1800"
# plot_dual_bound(organism, file_name)

organism = "iJR904"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name, csv=false)

# organism = "iJR904"
# file_name = "combinatorial_benders_1800"
# plot_dual_bound(organism, file_name)

organism = "iML1515"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name, csv=false)

# organism = "iML1515"
# file_name = "combinatorial_benders_1800"
# plot_dual_bound(organism, file_name)