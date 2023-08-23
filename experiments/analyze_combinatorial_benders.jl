using DataFrames
using CSV, JSON
using PyPlot

function plot_dual_bound(organism, file_name; fba_comparison=true)
    dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))

    objective_value = dict["objective_value"]
    objective_values = round.(dict["objective_values"], digits=4)
    solution = dict["solution"]
    time = dict["time"]
    termination = dict["termination"]
    iter = dict["iter"]

    # @show objective_value
    @show length(objective_values)
    @show objective_values[end]
    @show iter

    if fba_comparison
        # fba data 
        dict_ll_fba = JSON.parse(open("json/" * organism * "_loopless_fba_1800.json"))
    end

    fig = plt.figure(figsize=(6.5, 3.5))
    ax = fig.add_subplot(111)
    # ax.plot(collect(1:length(dual_bounds)), dual_bounds)
    ax.plot(collect(1:length(objective_values)), objective_values, label="combinatorial Benders'")
    # just the final objective value is available for loopless FBA, therefore it is mapped to the length
    # of the combinatorial Benders decomposition to compare the solutions 
    if fba_comparison
        ax.plot(collect(1:length(objective_values)), ones(length(objective_values))*round(dict_ll_fba["objective_value"],digits=5), label="ll-FBA")
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
    file = "plots/" * organism * "_" * file_name * ".pdf"
    savefig(file)
end 

# organism = "iAF692"
# file_name = "combinatorial_benders_fast_1800"
# plot_dual_bound(organism, file_name, fba_comparison=true)
# # plot_dual_bound(organism, file_name, fba_comparison=false)

# organism = "iJR904"
# file_name = "combinatorial_benders_fast_1800"
# plot_dual_bound(organism, file_name, fba_comparison=true)
# # plot_dual_bound(organism, file_name, fba_comparison=false)

# organism = "iML1515"
# file_name = "combinatorial_benders_fast_1800"
# plot_dual_bound(organism, file_name, fba_comparison=true)
# # plot_dual_bound(organism, file_name, fba_comparison=false)

organism = "Babjeviella_inositovora"
file_name = "combinatorial_benders_fast_7200"
plot_dual_bound(organism, file_name, fba_comparison=false)
