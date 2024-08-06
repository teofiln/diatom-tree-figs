library(treeio)
library(ggtree)
library(readxl)
library(tidyverse)
library(ape)
library(ggplot2)
library(wesanderson)

# read in tree file
diatom.tree <- treeio::read.mcmctree("data/summary-mean-brlens.nex")

# drop Bolidomonas
no_bolid <- treeio::drop.tip(
  object = diatom.tree,
  tip = diatom.tree@phylo$tip.label[1]
)

# read in metadata for tip labels
diatom.data <- read_excel("data/tree-figs-data.xlsx")

# get vector of tips to drop
##### --- prep for plotting --- #####
tips_only <- diatom.data[
  diatom.data$label %in% no_bolid@phylo$tip.label,
]
major_group_info <- split(tips_only$label, tips_only$clade.color)
major_groups <- groupOTU(no_bolid, major_group_info)

# rename tips
tips_only <- tips_only %>%
  unite("formatted.label", c(final.genus, final.species, strain),
        na.rm = TRUE,
        sep = " ")

# rename tips
major_groups@phylo <- rename_taxa(
  major_groups@phylo, tips_only,
  key = label, value = formatted.label
)

group_colors <- wesanderson::wes_palette(
  name = "Darjeeling1",
  n = 10,
  type = "continuous"
)
group_colors <- rev(group_colors)
group_colors <- c("darkgrey", "#5785c1", group_colors[-3])
group_colors[6] <- "#D37416"
group_colors[10] <- "#157764"

# uncomment to see colors
# "black"   "grey"    "#5785c1" "#5BBCD6" "#A1A376" "#F69100" "#F3A300" "#BCAA1E" "#50A45C" "#1C8E7A" "#8D473D" "#FF0000"

time.labels <- as.character(seq(300, 0, by=-50))

p1 <-
  ggtree(
    major_groups,
    size = 0.7, aes(color = group)
    ) +
  geom_range("height_0.95_HPD", color = "#4361e7", size = 1.5, alpha = .7) +
  geom_tiplab(size = 2.5, family = "Helvetica") +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_fill_manual(values = group_colors) +
  theme_tree2() +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.35, vjust = -0.1)) +
  xlab("Million years ago (Ma)") +
  theme(axis.line.x=element_line(), axis.line.y=element_blank()) +
  theme(panel.grid.major.x=element_line(color="darkgrey", linetype="dashed", linewidth = .4),
        panel.grid.major.y=element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  scale_y_continuous(limits = c(0, 258), expand = expansion(add = c(2, 0)))
  
p2 <- revts(p1)
p3 <- p2 + scale_x_continuous(limits = c(-300, 130), breaks = seq(-300, 0, by=50), labels = time.labels)

ggsave(
  p3,
  file = "full-tree-figure/full-tree.pdf", width = 8, height = 24
)
 
