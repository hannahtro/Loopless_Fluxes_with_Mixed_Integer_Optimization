import numpy as np               
import matplotlib.pyplot as plt    
import matplotlib.patches as patches

# Select length of axes and the space between tick labels
xmin, xmax, ymin, ymax = 0, 3.5, 0, 3.5
ticks_frequency = 1

# Plot points
fig, ax = plt.subplots(figsize=(10, 10))
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
ax.set_xlabel('x', size=14, labelpad=-24, x=1.03)
ax.set_ylabel('y', size=14, labelpad=-21, y=1.02, rotation=0)

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

# feasible region
coord = [[2,3.5], [3,3], [3,0.5], [0,0.5]]
coord.append(coord[0])
xs, ys = zip(*coord)
ax.plot(xs, ys, color='cornflowerblue', label='constraints') 

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

# integer points 
plt.plot(1, 1, 'ro', color='cornflowerblue') 
plt.plot(1, 2, 'ro', color='cornflowerblue') 
plt.plot(2, 1, 'ro', color='cornflowerblue') 
plt.plot(3, 1, 'ro', color='cornflowerblue') 
plt.plot(2, 2, 'ro', color='cornflowerblue') 
plt.plot(3, 2, 'ro', color='cornflowerblue') 
plt.plot(2, 3, 'ro', color='cornflowerblue') 
plt.plot(3, 3, 'ro', color='cornflowerblue') 

# S-free set
circle = plt.Circle((2, 3.5), 0.5, color='r', fill=False)
# rect = patches.Rectangle((-1, 0), 7, 6, linewidth=1, edgecolor='orange', facecolor='orange', alpha=0.3, label='S-free set')
ax.add_patch(circle)

# plot conic relaxation
x1, y1 = [2, 5], [3.5,2]
x2, y2 = [2, -0.33], [3.5, 0]
ax.plot(x1, y1, color='r', linestyle='dotted', alpha=1, label='conic relaxation')
ax.plot(x2, y2, color='r', linestyle='dotted', alpha=1)

# intersection cut 
x3, y3 = [1.723, 2.447], [3.084, 3.278]
ax.plot(x3, y3, color='g', linestyle='-', alpha=1, label='intersection cut')
plt.plot(1.723, 3.084, 'ro', color='g') 
plt.plot(2.447, 3.278, 'ro', color='g') 
line(0.268, 2.62, color="g")

# optimal point
plt.plot(2, 3.5, 'ro', color='r', alpha=1, label='optimal solution') 

# optimal point after cut
# plt.plot(2.45, 3.275, 'ro', color='r', alpha=1) 

fig.legend(loc='center', bbox_to_anchor=(0.52, 0.05), ncol=2)

plt.savefig("plots/intersection_cut.pdf", format="pdf", bbox_inches="tight")