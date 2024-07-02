# Fan tree with time scale

library(ggplot2)
library(ggtree)
library(treeio)
library(deeptime)
library(ape)
library(readxl)
library(viridis)
library(dplyr)
library(purrr)
library(wesanderson)

##### --- the tree --- #####
diatom.tree <- treeio::read.mcmctree("fan-figure/summary-topology.nex")
no_bolid <- treeio::drop.tip(
  object = diatom.tree,
  tip = diatom.tree@phylo$tip.label[1]
)

##### --- the data --- #####
diatom.data <- readxl::read_excel("fan-figure/voucher-list.xlsx")

##### --- prep for plotting --- #####
tips_only <- diatom.data[
  diatom.data$label %in% no_bolid@phylo$tip.label,
]
major_group_info <- split(tips_only$label, tips_only$clade.color)
major_groups <- groupOTU(no_bolid, major_group_info)

##### --- data frame for time scale labels --- #####
periods <- data.frame(
  start = c(0, 23, 66, 145, 201, 252, 298),
  end = c(23, 66, 145, 201, 252, 298, 359),
  mid = c(11.5, 44.5, 105.5, 173, 226.5, 275, 328.5),
  period = c(
    "Neogene",
    "Paleogene",
    "Cretaceous",
    "Jurassic",
    "Triassic",
    "Permian",
    "Carboniferous"
  ),
  abbrv = c("Ng", "Pg", "K", "J", "Tr", "P", "C")
)

periods <- periods[-7, ]

##### --- data frame for clade labels --- #####
no_bolid_phy <- no_bolid@phylo
no_bolid_dat <- diatom.data[diatom.data$label %in% no_bolid_phy$tip.label, ]


fan_timescale_plot <- function(wes_pal = "Darjeeling1") {
  n_cols <- length(unique(no_bolid_dat$clade.color))
  bar_colors <- wesanderson::wes_palette(
    name = wes_pal,
    n = n_cols,
    type = "continuous"
  )
  bar_colors <- rev(bar_colors)

  text_colors <- rev(c(rep("#FFFFFF", 4), rep("#000000", 6)))

  clade_mrcas <- no_bolid_dat |>
    dplyr::group_split(clade.color) |>
    purrr::map_dfr(~ {
      anc <- ape::getMRCA(no_bolid_phy, .x$label)
      data.frame(group = unique(.x$clade.color), node = anc)
    }) |>
    dplyr::mutate(label = factor(dplyr::row_number())) |>
    dplyr::mutate(clade_color = bar_colors)

  ##### --- plot --- #####
  p <- revts(ggtree(major_groups, aes(color = group))) +
    coord_geo_radial(
      dat = "periods",
      start = -0.5 * pi,
      end = 0.5 * pi,
      fill = c("#FFFFFF", "#f7f7f7"),
      color = "darkgrey",
      lty = "dashed",
      lwd = 0.15,
    ) +
    geom_range("height_0.95_HPD", color = "#4361e7", size = 0.75, alpha = .7) +
    ggtree::geom_cladelab(
      data = clade_mrcas,
      mapping = aes(
        node = node, label = label, fill = label,
        hjust = 0.5, yjust = 1
      ),
      geom = "label",
      show.legend = FALSE,
      extend = 0.25,
      offset = 7,
      offset.text = 15,
      barsize = 1.5,
      barcolour = bar_colors,
      textcolour = text_colors,
      label.r = unit(0.5, "lines")
    ) +
    ggplot2::geom_text(
      data = periods,
      ggplot2::aes(x = -mid, y = 257, label = abbrv),
      inherit.aes = FALSE,
      size = 3,
      color = "#797979"
    ) +
    ggplot2::scale_color_manual(values = c("black", bar_colors)) +
    ggplot2::scale_fill_manual(values = c(bar_colors)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(add = c(6, 0)),
      guide = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(-300, 0, 50),
      labels = abs(seq(-300, 0, 50)),
      expand = ggplot2::expansion(add = c(NA, 0)),
    ) +
    ggplot2::labs(y = "Time (Ma)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(vjust = -2, hjust = 0.25, size = 10),
    )

  file_pdf <- sprintf("fan-figure/fan-timescale-%s.png", wes_pal)
  file_png <- sprintf("fan-figure/fan-timescale-%s.pdf", wes_pal)

  ggplot2::ggsave(file_pdf, p, width = 10, height = 7, dpi = 300)
  ggplot2::ggsave(file_png, p, width = 10, height = 7)

  # needs pkg `svglite`
  # ggplot2::ggsave("fan-timescale.svg", p, width = 10, height = 7)
}

all_wes_palettes <- names(wesanderson::wes_palettes)

purrr::map(all_wes_palettes, ~ fan_timescale_plot(.x))
