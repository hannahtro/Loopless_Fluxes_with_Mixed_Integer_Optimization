import numpy as np               
import matplotlib.pyplot as plt    
import matplotlib.patches as patches

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "lmodern"
})
plt.rcParams.update({'figure.autolayout': True})

# Select length of axes and the space between tick labels
xmin, xmax, ymin, ymax = -3, 6, -2, 5
ticks_frequency = 1

# Plot points
fig, ax = plt.subplots(figsize=(12, 6))
# ax.scatter(xs, ys, c=colors)

# Set identical scales for both axes
ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin-1, ymax+1), aspect='equal')

# Set bottom and left spines as x and y axes of coordinate system
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Create 'x' and 'y' labels placed at the end of the axes
ax.set_xlabel(r"$\Delta \mu_i h$", size=14, labelpad=-24, x=1.05)
ax.set_ylabel(r"$v_i$", size=14, labelpad=-21, y=1.02, rotation=0)

# Create custom major ticks to determine position of tick labels
x_ticks = np.arange(xmin, xmax+1, ticks_frequency)
y_ticks = np.arange(ymin, ymax+1, ticks_frequency)
ax.set_xticks(x_ticks[x_ticks != 0])
ax.set_yticks(y_ticks[y_ticks != 0])

# Create minor ticks placed at each integer to enable drawing of minor grid
# lines: note that this has no effect in this example with ticks_frequency=1
ax.set_xticks(np.arange(xmin, xmax+1), minor=True)
ax.set_yticks(np.arange(ymin, ymax+1), minor=True)

# Draw major and minor grid lines
ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

# Draw arrows
arrow_fmt = dict(markersize=4, color='black', clip_on=False)
ax.plot((1), (0), marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
ax.plot((0), (1), marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

# set S
rect = patches.Rectangle((-6, 0), 5, 6, linewidth=1, edgecolor='black', facecolor='b', alpha=0.2, label='thermodynamical feasible region')
ax.add_patch(rect)
rect = patches.Rectangle((1, -6), 6, 6, linewidth=1, edgecolor='black', facecolor='b', alpha=0.2)
ax.add_patch(rect)

# feasible region 
plt.rcParams.update({'hatch.color': 'gray'})
coord = [[1.0,0], [3.5,0], [3,-2], [1.0,-2]]
coord.append(coord[0])
xs, ys = zip(*coord)
ax.plot(xs, ys, color='lightgray') 
ax.fill(xs, ys, "lightgray", alpha=0.5, hatch="/")

coord = [[-1.0,0], [-2,0], [-2,4], [-1.0,4]]
coord.append(coord[0])
xs, ys = zip(*coord)
ax.plot(xs, ys, color='lightgray') 
ax.fill(xs, ys, "lightgray", alpha=0.5, hatch="/")

# set P
coord = [[-2,4], [2,4], [4,2], [3,-2], [-2,-2]]
coord.append(coord[0])
xs, ys = zip(*coord)
ax.plot(xs, ys, color='cornflowerblue', label='constraints') 

# S-free set
rect = patches.Rectangle((-1, 0), 8, 6, linewidth=1, edgecolor='orange', facecolor='orange', alpha=0.3, label='S-free set')
ax.add_patch(rect)

# plot conic relaxation
x1, y1 = [2, 7], [4, -1]
x2, y2 = [2, -6], [4, 4]
ax.plot(x1, y1, color='g', linestyle='dotted', alpha=1, label='conic relaxation')
ax.plot(x2, y2, color='g', linestyle='dotted', alpha=1)

# intersection cut 
x3, y3 = [-4.5, 7], [6, -0.6]
ax.plot(x3, y3, color='r', linestyle='-', alpha=1, label='intersection cut')

# optimal point
plt.plot(2, 4, 'ro', color='r', alpha=1, label='optimal solution') 

# fig.legend(loc='center', bbox_to_anchor=(0.52, 0.05), ncol=3)
plt.savefig("plots/intersection_cut.pdf", format="pdf", bbox_inches="tight")