import numpy as np               
import matplotlib.pyplot as plt    
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator

dp = False
big_m = False
convex_hull = True 

plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern"
    })

# Select length of axes and the space between tick labels
xmin, xmax, ymin, ymax = 0.5, 3, 0.5, 3
ticks_frequency = 1

# Plot points
fig, ax = plt.subplots(figsize=(5, 5))
# ax.scatter(xs, ys, c=colors)

# Set identical scales for both axes
# ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin-1, ymax+1), aspect='equal')
ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin-1, ymax+1))

# Set bottom and left spines as x and y axes of coordinate system
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Create 'x' and 'y' labels placed at the end of the axes
ax.set_xlabel('x', size=12, labelpad=-24, x=1.03)
ax.set_ylabel('y', size=12, labelpad=-21, y=1.02, rotation=0)

# Create custom major ticks to determine position of tick labels
x_ticks = np.arange(0, 4, ticks_frequency)
y_ticks = np.arange(0, 4, ticks_frequency)
ax.set_xticks(x_ticks[x_ticks != 0])
ax.set_yticks(y_ticks[y_ticks != 0])
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# Create minor ticks placed at each integer to enable drawing of minor grid
# lines: note that this has no effect in this example with ticks_frequency=1
# ax.set_xticks(np.arange(xmin, xmax+1), minor=True)
# ax.set_yticks(np.arange(ymin, ymax+1), minor=True)

# Draw major and minor grid lines
ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

# Draw arrows
arrow_fmt = dict(markersize=4, color='black', clip_on=False)
ax.plot((1), (0), marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
ax.plot((0), (1), marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

# feasible region
coord = [[1.0,2.0], [1.0,0.5], [0,0.5], [1.0,2.0]]
coord.append(coord[0])
xs, ys = zip(*coord)
ax.plot(xs, ys, color='cornflowerblue', label='constraints') 
ax.fill(xs, ys, "cornflowerblue", alpha=0.8)

coord = [[2,3.5], [3,3], [3,2.0], [2,2.0]]
coord.append(coord[0])
xs, ys = zip(*coord)
ax.plot(xs, ys, color='orange', label='constraints')

if dp:
    ax.fill(xs, ys, "orange", alpha=0.5)
else:
    ax.fill(xs, ys, "orange", alpha=0.8)
        
if convex_hull:
    # convex hull 
    plt.rcParams.update({'hatch.color': 'gray'})
    coord = [[1.0,0.5], [0,0.5], [2,3.5], [3,3], [3,2.0]]
    coord.append(coord[0])
    xs, ys = zip(*coord)
    ax.plot(xs, ys, color='lightgray', label='convex hull') 
    ax.fill(xs, ys, "lightgray", alpha=0.5, hatch="+")

if big_m:
    # big-M 
    plt.rcParams.update({'hatch.color': 'gray'})
    coord = [[1.0,0.5], [0,0.5], [2,3.5], [3,3], [3,2.0], [3,0.5]]
    coord.append(coord[0])
    xs, ys = zip(*coord)
    # ax.plot(xs, ys, color='lightgray', label='convex hull') 
    ax.fill(xs, ys, "lightgray", alpha=0.5, hatch="+")

def line(slope, intercept, color="cornflowerblue"):
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--', color=color)

def y_line(bound, color='cornflowerblue'):
    axes = plt.gca()
    y_vals = np.array(axes.get_ylim())
    x_vals = [bound for val in y_vals]
    plt.plot(x_vals, y_vals, '--', color=color)

# line(0, 3.5)
# line(0, 0.5)
# line(1.5, 0.5)
# line(-0.5, 4.5)
# y_line(3)
# y_line(0)

# # integer points 
# plt.plot(1, 1, 'ro', color='cornflowerblue') 
# plt.plot(1, 2, 'ro', color='cornflowerblue') 
# plt.plot(2, 1, 'ro', color='cornflowerblue') 
# plt.plot(3, 1, 'ro', color='cornflowerblue') 
# plt.plot(2, 2, 'ro', color='cornflowerblue') 
# plt.plot(3, 2, 'ro', color='cornflowerblue') 
# plt.plot(2, 3, 'ro', color='cornflowerblue') 
# plt.plot(3, 3, 'ro', color='cornflowerblue') 

# optimal point
plt.plot(2, 3.5, 'ro', color='r', alpha=1, label='optimal solution') 

# optimal point after cut
# plt.plot(2.45, 3.275, 'ro', color='r', alpha=1) 

# fig.legend(loc='center', bbox_to_anchor=(0.52, 0.05), ncol=2)

if dp:
    plt.savefig("plots/dp.pdf", format="pdf", bbox_inches="tight")
elif convex_hull:
    plt.savefig("plots/dp_convex_hull.pdf", format="pdf", bbox_inches="tight")
elif big_m:
    plt.savefig("plots/dp_big_m.pdf", format="pdf", bbox_inches="tight")