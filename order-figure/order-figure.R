library(treeio)
library(ggtree)
library(readxl)
library(dplyr)
library(ape)
library(ggplot2)
library(wesanderson)
library(here)

tree.file     <- here("data", "thal-raphid-tax0.6-pmsf.treefile")
metadata.file <-  here("data", "tree-figs-data.xlsx")
metadata.tree <- here("data", "tree-figs-data.xlsx")

tree <- ape::read.tree(tree.file)

# read in metadata for tip labels
metadata <- read_excel(metadata.file)

# read in metadata for making tree figs -- this will allow us to count orders that were
# plotted twice, meaning they are not monophyletic
tree.figs.metadata <- read_excel(metadata.tree)

# root the tree
rooted.tree <- root.phylo(
  phy = tree,
  outgroup = "Aureococcusanoph",
  resolve.root = TRUE,
  edgelabel = TRUE
)

# get vector of tips to drop
drop.dat <- metadata |> dplyr::filter(order.tree == "drop")
drop.vec <- drop.dat$label

# drop tips
rooted.tree <- drop.tip(rooted.tree, tip = drop.vec)

# copy metadata object to metadata_tips_only with only the rows matching
# rooted.tree$tip.label, then make a list of lists from that to input to groupOTU
metadata_tips_only <- metadata[metadata$label %in% rooted.tree$tip.label, ]

### group by morphological "group" == radial, polar, araphid, or raphid
# creates a list of lists, each one is the tips that belong to each group
major_group_info <- split(
  metadata_tips_only$label, metadata_tips_only$clade.color
)
major_groups <- groupOTU(rooted.tree, major_group_info)

# rename tips
major_groups.tree <- rename_taxa(
  major_groups, metadata_tips_only,
  key = label, value = order
)

group_colors <- wesanderson::wes_palette(
  name = "Darjeeling1",
  n = 10,
  type = "continuous"
)
group_colors <- rev(group_colors)
group_colors <- c("black", "grey", "#5785c1", group_colors[-3])

# uncomment to see colors
# "black"   "grey"    "#5785c1" "#5BBCD6" "#A1A376" "#F69100" "#F3A300" "#BCAA1E" "#50A45C" "#1C8E7A" "#8D473D" "#FF0000"

# manually modify a few of the colors to increase contrast of adjacent colors
group_colors[6] <- "#D37416"
group_colors[10] <- "#157764"


major_group_tree_fig <-
  ggtree(
    major_groups.tree,
    size = 0.4, aes(color = group), branch.length = "none"
  ) +
  ggplot2::xlim(0, 40) +
  geom_tiplab(size = 2, family = "Helvetica") +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_fill_manual(values = group_colors)

outfile = here("order-figure", "orders.pdf")

ggsave(
  major_group_tree_fig,
  file = outfile, width = 2.75, height = 8.73
)
