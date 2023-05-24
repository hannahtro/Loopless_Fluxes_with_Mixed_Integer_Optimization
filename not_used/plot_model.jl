using Escher, CairoMakie, ColorSchemes

f = Figure(resolution = (1200, 800));
ax = Axis(f[1, 1]);
###### PLOT FUNCTION
hidexdecorations!(ax)
hideydecorations!(ax)
escherplot!(ax, "data/e_coli_core_map.json")

save("data/e_coli_core_map.png",f)