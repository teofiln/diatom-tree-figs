library(here)
library(readxl)
library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
library(phytools)
library(phangorn)
library(TreeTools)
library(deeptime)

# read in species tree with divergence times
diatom.tree <- read.mcmctree('summary-topology.nex')

# read in metadata for tip labels
metadata <- read_excel(here('cultures-metadata', 'voucher-list.xlsx'), sheet = '292-dataset')

# make a column with formatted tip label
metadata <- metadata %>%
  unite("formatted.label", c(final.genus, final.species),
        na.rm = TRUE,
        sep = " ")

# reduce metadata to include only rows that correspond to species present in the tree
# not doing this will mess up groupOTU
metadata_tips_only <- metadata[metadata$label %in% diatom.tree@phylo$tip.label,]

# group by clade (need to make a column for this in metadata)
# major_group_info <- split(metadata_tips_only$label, metadata_tips_only$clade)
# major_groups <- groupOTU(diatom.tree, major_group_info)

# rename tips
# major_groups.tree <- rename_taxa(major_groups, metadata_tips_only, key = label, value = formatted.label)

tree.fig <- ggtree(major_groups.tree, size = 0.5, 
                               layout="fan", open.angle = 180)
                               # aes(color=group))

ggsave(tree.fig, file='test.pdf', width = 11, height = 8)


########## CODE FROM BIRD PAPER ##########
# Half-circular tree
## To add custom annotations of the circular time scale make a df
time <- data.frame(max_age = c(66, 23.03, 20, 40, 60, 80, 100),
                   min_age = c(66, 23.03, 0,  20, 40, 60, 80))
# tree <- groupOTU(tree, cls)
ggtree(tree, size=1, aes(color=group), layout="fan", open.angle = 180) + theme_tree2() + 
  # add nodes for fossils
  geom_nodepoint(aes(shape=fossil), size=0.5, color="grey20") +
  scale_shape_manual(values=c(NA, 16)) +
  theme(legend.position = "none") +
  coord_geo_polar(dat = "periods", fill=c("white"), color="grey80", lty=1) +
  # Show time slots every 20 my and also K-Pg and Pg-Ng boundaries
  # add CIs
  geom_range(range="height_0.95_HPD", color="grey70", alpha=0.7, size=0.5) +
  theme(plot.margin = unit(c(1,1,0,-5), "mm"), legend.position="none",
        axis.text.r = element_text(size = 2)) +
  scale_x_continuous(breaks = seq(-100, 0, 20),
                     labels = seq(100, 0, -20)) +
  scale_color_manual(values=c("grey20", unique(a$color_mag7_hex))) -> p3
revts(p3) -> p3
