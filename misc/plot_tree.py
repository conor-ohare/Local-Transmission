import dendropy
from Bio import Phylo
import matplotlib.pyplot as plt

# Load the tree using DendroPy
tree_file = "data/NEA_D4.tree"
tree = dendropy.Tree.get(path=tree_file, schema="nexus")

# Save the tree in Newick format to read it with Bio.Phylo
tree.write(path="NEA_D4.newick", schema="newick")

# Read the Newick file using Bio.Phylo
phylo_tree = Phylo.read("NEA_D4.newick", "newick")

# Plot the tree with branch lengths
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(phylo_tree, do_show=False, branch_labels=lambda c: f"{c.branch_length:.2f}" if c.branch_length else "")

# Save the plot
plt.savefig("NEA_D4_tree_plot.png", dpi=300)
print("Tree plot saved as 'NEA_D4_tree_plot.png'")

# Display the plot
plt.show()