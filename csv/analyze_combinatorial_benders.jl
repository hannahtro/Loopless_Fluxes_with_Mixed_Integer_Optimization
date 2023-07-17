using DataFrames
using CSV
using PyPlot

function plot_dual_bound(organism, file_name; fba_comparison=true)
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
    @show objective_values
    @show iter

    if fba_comparison
        # fba data 
        df_fba = first(CSV.read(organism * "/server/" * organism * "_loopless_fba_1800.csv", DataFrame),1)
    end

    fig = plt.figure(figsize=(6.5, 3.5))
    ax = fig.add_subplot(111)
    # ax.plot(collect(1:length(dual_bounds)), dual_bounds)
    ax.plot(collect(1:length(objective_values)), objective_values, label="combinatorial Benders'")
    # just the final objective value is available for loopless FBA, therefore it is mapped to the length
    # of the combinatorial Benders decomposition to compare the solutions 
    if fba_comparison
        ax.plot(collect(1:length(objective_values)), ones(length(objective_values))*df_fba[!,:objective_value][1], label="FBA")
        # ax.plot(collect(1:length(df_fba[!,:objective_value])), df_fba[!,:objective_value], label="FBA")
    end

    ylabel("dual bound")
    #locator_params(axis="y", nbins=4)
    xlabel("iteration")
    ax.grid()
    ax.legend()
    fig.tight_layout()
    f = matplotlib.ticker.FormatStrFormatter("%1.1f") # Define format of tick labels
    ax.xaxis.set_major_formatter(f) # Set format of tick labels

    if fba_comparison
        file_name = file_name * "_fba_comparison"
    end
    file = "../plots/" * organism * "_" * file_name * ".pdf"
    savefig(file)
end 

organism = "iAF692"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name, fba_comparison=true)
# plot_dual_bound(organism, file_name, fba_comparison=false)

# organism = "iAF692"
# file_name = "combinatorial_benders_1800"
# plot_dual_bound(organism, file_name)

organism = "iJR904"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name, fba_comparison=true)
# plot_dual_bound(organism, file_name, fba_comparison=false)

# organism = "iJR904"
# file_name = "combinatorial_benders_1800"
# plot_dual_bound(organism, file_name)

organism = "iML1515"
file_name = "combinatorial_benders_fast_1800"
plot_dual_bound(organism, file_name, fba_comparison=true)
# plot_dual_bound(organism, file_name, fba_comparison=false)

# organism = "iML1515"
# file_name = "combinatorial_benders_1800"
# plot_dual_bound(organism, file_name)