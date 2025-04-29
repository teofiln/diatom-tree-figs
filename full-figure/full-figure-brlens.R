library(treeio)
library(ggtree)
library(readxl)
library(tidyverse)
library(ape)
library(ggplot2)
library(wesanderson)
library(here)

# read in tree file
diatom.tree <- treeio::read.tree("data/thal-raphid-tax0.6-pmsf.treefile")

# read in metadata for tip labels
diatom.data <- read_excel("data/tree-figs-data.xlsx")

# root the tree
rooted.tree <- root.phylo(
  phy = diatom.tree,
  outgroup = "Aureococcusanoph",
  resolve.root = TRUE,
  edgelabel = TRUE
)

# get vector of tips to drop
drop.dat <- diatom.data |> dplyr::filter(keep == "drop")
drop.vec <- drop.dat$label

# drop tips
rooted.tree <- drop.tip(rooted.tree, tip = drop.vec)

# drop non-Nagoya taxa and Aureococcus
rooted.tree <- drop.tip(
  rooted.tree,
  tip = c("Naviculoid_ECT2AJA-163_L210", "Parlibellus_cf_delognei__ECT2AJA-152_L165", "Aureococcusanoph")
)

##### --- prep for plotting --- #####
tips_only <- diatom.data[
  diatom.data$label %in% rooted.tree$tip.label,
]
major_group_info <- split(tips_only$label, tips_only$clade.color)
major_groups <- groupOTU(rooted.tree, major_group_info)

# create formatted tip labels
tips_only <- tips_only %>%
  unite("formatted.label", c(final.genus, final.species, strain),
        na.rm = TRUE,
        sep = " ")

# rename tips
major_groups <- rename_taxa(
  major_groups, tips_only,
  key = label, value = formatted.label
)

# define color palette (borrowed from order-figure.R)
group_colors <- wesanderson::wes_palette(
  name = "Darjeeling1",
  n = 10,
  type = "continuous"
)
group_colors <- rev(group_colors)
group_colors <- c("darkgrey", "darkgrey", "#5785c1", group_colors[-3])
group_colors[6] <- "#D37416"
group_colors[10] <- "#157764"

# # define color palette
# group_colors <- wesanderson::wes_palette(
#   name = "Darjeeling1",
#   n = 10,
#   type = "continuous"
# )
# group_colors <- rev(group_colors)
# group_colors <- c("grey", "grey", "#5785c1", group_colors[-3])
# group_colors[6] <- "#D37416"
# group_colors[10] <- "#157764"

# # define color palette
# group_colors <- wesanderson::wes_palette(
#   name = "Darjeeling1",
#   n = 11,
#   type = "continuous"
# )
# group_colors <- rev(group_colors)
# group_colors <- c("darkgrey", "#5785c1", group_colors[-3])
# group_colors[6] <- "#D37416"
# group_colors[10] <- "#157764"

##### --- make tree plot --- #####
p1 <-
  ggtree(
    major_groups,
    size = 0.7, aes(color = group)
    ) +
  geom_tiplab(size = 2.5, family = "Helvetica") +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_fill_manual(values = group_colors) +
  geom_treescale(x=0, y=40, width=0.1, color='black', linesize = 1, fontsize = 6, offset = 1.5) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  ggplot2::xlim(0, 3.5)

#  theme_tree2() +
#  theme(axis.line.x=element_line(), axis.line.y=element_blank()) +
#  scale_y_continuous(expand = expansion(add = c(2, 1)))
  
ggsave(
  p1,
  file = "full-figure/full-tree-brlens.R2.pdf", width = 8, height = 25
)

