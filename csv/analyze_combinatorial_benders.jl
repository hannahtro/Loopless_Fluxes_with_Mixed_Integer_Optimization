using DataFrames
using CSV
using PyPlot

function plot_dual_bound(organism, file_name)
    df = first(CSV.read(organism * "_" * file_name * ".csv", DataFrame),1)
    @show df

    objective_value = df[!,:objective_value]
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
    @show dual_bounds[end]
    @show length(dual_bounds)
    @show iter

    fig = plt.figure(figsize=(6.5, 3.5))
    ax = fig.add_subplot(111)
    ax.plot(collect(1:length(dual_bounds)), dual_bounds)

    ylabel("dual bound")
    #locator_params(axis="y", nbins=4)
    xlabel("iteration")
    ax.grid()
    fig.tight_layout()
    f = matplotlib.ticker.FormatStrFormatter("%1.1f") # Define format of tick labels
    ax.xaxis.set_major_formatter(f) # Set format of tick labels

    file = "../plots/" * organism * "_" * file_name * ".pdf"
    savefig(file)
end 

organism = "iAF692"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name)

organism = "iAF692"
file_name = "combinatorial_benders_1800"
plot_dual_bound(organism, file_name)

organism = "iJR904"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name)