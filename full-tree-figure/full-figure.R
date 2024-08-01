library(treeio)
library(ggtree)
library(readxl)
library(tidyverse)
library(ape)
library(ggplot2)
library(wesanderson)

# read in tree file
diatom.tree <- treeio::read.mcmctree("data/summary-mean-brlens.nex")

# read in metadata for tip labels
metadata <- read_excel("data/tree-figs-data.xlsx")

# root the tree
rooted.tree <- root.phylo(
  phy = diatom.tree@phylo,
  outgroup = "Bolidomonas_pacifica_CCMP1866",
  resolve.root = TRUE,
  edgelabel = TRUE
)

# get vector of tips to drop
drop.dat <- metadata |> dplyr::filter(keep == "drop")
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
metadata_tips_only <- metadata_tips_only %>%
  unite("formatted.label", c(final.genus, final.species, strain),
        na.rm = TRUE,
        sep = " ")

# rename tips
major_groups.tree <- rename_taxa(
  major_groups, metadata_tips_only,
  key = label, value = formatted.label
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
    size = 0.4, aes(color = group)
    ) +
  scale_y_reverse() +
#  ggplot2::xlim(0, 40) +
  geom_tiplab(size = 2, family = "Helvetica") +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_fill_manual(values = group_colors)

outfile = here("full-tree-figure/full-tree.pdf")

ggsave(
  major_group_tree_fig,
  file = outfile, width = 6, height = 20
)

