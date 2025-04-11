import matplotlib.pyplot as plt
# Adjusted version with arrows placed between nodes, not inside
fig, ax = plt.subplots(figsize=(5, 2))
ax.set_xlim(-0.5, 2.5)
ax.set_ylim(-0.5, 0.5)
ax.axis("off")

# Node positions
nodes = {"S": (0, 0), "I": (1, 0), "R": (2, 0)}
colors = {"S": "deepskyblue", "I": "red", "R": "lightgreen"}

# Draw nodes
for node, (x, y) in nodes.items():
    ax.add_patch(plt.Circle((x, y), 0.2, color=colors[node], ec="black", lw=2))
    ax.text(x, y, node, ha="center", va="center", fontsize=16, fontweight="bold", color="black")

# Draw arrows between nodes, placed outside of the circles
arrow_params = dict(arrowstyle="->", color="black", lw=2)
ax.annotate("", xy=(0.7, 0), xytext=(0.35, 0), arrowprops=arrow_params)  # S -> I
ax.annotate("", xy=(1.7, 0), xytext=(1.35, 0), arrowprops=arrow_params)  # I -> R

# Death arrows, placed just below each node
for x, y in nodes.values():
    ax.annotate("", xy=(x, y - 0.5), xytext=(x, y - 0.25), arrowprops=arrow_params)  # Downward death arrows

# Birth arrow, placed before S
ax.annotate("", xy=(-0.3, 0), xytext=(-0.6, 0), arrowprops=arrow_params)  # Birth into S

plt.show()
