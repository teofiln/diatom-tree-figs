library(treeio)
library(ggtree)
library(readxl)
library(tidyverse)
library(ape)
library(ggplot2)
library(wesanderson)
library(here)

# read in tree file
diatom.tree <- treeio::read.mcmctree("data/summary-mean-brlens.nex")

# drop Bolidomonas
no_bolid <- treeio::drop.tip(
  object = diatom.tree,
  tip = diatom.tree@phylo$tip.label[1]
)

# read in metadata for tip labels
diatom.data <- read_excel("data/tree-figs-data.xlsx")

##### --- prep for plotting --- #####
tips_only <- diatom.data[
  diatom.data$label %in% no_bolid@phylo$tip.label,
]
major_group_info <- split(tips_only$label, tips_only$clade.color)
major_groups <- groupOTU(no_bolid, major_group_info)

n1 <- getMRCA(phy=major_groups@phylo, 
              c('Lithodesmium_intricatum_ECT2AJA-029_L217', 'Discostella_pseudostelligera_AJA075-4_R26'))
n2 <- getMRCA(phy=major_groups@phylo,
              c('Toxarium_undulatum_ECT3802_L241', 'Orthoseira_roeseana_CBG002_L503'))
n3 <- getMRCA(phy=major_groups@phylo,
              c('Toxarium_undulatum_ECT3802_L241', 'Lampriscus_shadboltianum_ECT2AJA-054_R71'))
n4 <- getMRCA(phy=major_groups@phylo,
              c('Attheya_longicornis_ECT2AJA-053_R15', 'Terpsinoe_americana_ECT2AJA-024_L208'))
n5 <- getMRCA(phy=major_groups@phylo,
              c('Toxarium_undulatum_ECT3802_L241', 'Papiliocellulus_elegans_CCMP3125_L103'))
n6 <- getMRCA(phy=major_groups@phylo,
              c('Lithodesmium_intricatum_ECT2AJA-029_L217', 'Toxarium_undulatum_ECT3802_L241'))
n7 <- getMRCA(phy=major_groups@phylo,
              c('Attheya_longicornis_ECT2AJA-053_R15', 'Eunotia_naegelii_FD354_L122'))
# n8 <- getMRCA(phy=major_groups@phylo,
#               c())

# create formatted tip labels
tips_only <- tips_only %>%
  unite("formatted.label", c(final.genus, final.species, strain),
        na.rm = TRUE,
        sep = " ")

# rename tips
major_groups@phylo <- rename_taxa(
  major_groups@phylo, tips_only,
  key = label, value = formatted.label
)

# define color palette
group_colors <- wesanderson::wes_palette(
  name = "Darjeeling1",
  n = 10,
  type = "continuous"
)
group_colors <- rev(group_colors)
group_colors <- c("darkgrey", "#5785c1", group_colors[-3])
group_colors[6] <- "#D37416"
group_colors[10] <- "#157764"

# labels for timescale
time.labels <- as.character(seq(300, 0, by=-25))

##### --- make tree plot --- #####
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
  theme(panel.grid.major.x=element_line(color="darkgrey", linetype="dashed", linewidth = .3),
        panel.grid.major.y=element_blank()) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  scale_y_continuous(expand = expansion(add = c(2, 1)))
  
p2 <- revts(p1)
p3 <- p2 + scale_x_continuous(limits = c(-300, 100), breaks = seq(-300, 0, by=25), labels = time.labels)
# thal + litho clade
p4 <- p3 + geom_point2(aes(subset=node==n1), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-9))
# chaetocerotales + sister clade
p4 <- p4 + geom_point2(aes(subset=node==n2), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-7))
# toxariales clade
p4 <- p4 + geom_point2(aes(subset=node==n3), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-9))
# attheya clade
p4 <- p4 + geom_point2(aes(subset=node==n4), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-15))
# toxariales + sister clade
p4 <- p4 + geom_point2(aes(subset=node==n5), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-9))
# mediophytes minus attheya clade
p4 <- p4 + geom_point2(aes(subset=node==n6), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-8))
# attheya plus pennates
p4 <- p4 + geom_point2(aes(subset=node==n7), 
                       color='black', size=3, shape = 19, position = position_nudge(x=-9))


ggsave(
  p4,
  file = "full-figure/full-tree.pdf", width = 8, height = 25
)

