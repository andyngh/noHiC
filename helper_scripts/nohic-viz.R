#!/usr/bin/env Rscript

ensure_pkgs <- function(pkgs) {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (user_lib == "") {
    user_lib <- file.path(path.expand("~"), "R", "library")
    Sys.setenv(R_LIBS_USER = user_lib)
  }
  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org", lib = user_lib, dependencies = TRUE)
  }
}
ensure_pkgs(c("ggplot2", "readr", "dplyr"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: nohic-viz.R <chrom_lengths.csv> <misassemblies.csv> <output_prefix> <scale_factor>")
}
chrom_csv <- args[1]
mis_csv   <- args[2]
out_pref  <- args[3]
scale_factor <- args[4]

chroms <- readr::read_csv(chrom_csv, show_col_types = FALSE, col_names = TRUE)
if (!all(c("chrom", "length") %in% names(chroms))) {
  if (ncol(chroms) >= 2) {
    names(chroms)[1:2] <- c("chrom","length")
  } else stop("First CSV must have two columns: chrom,length (with header).")
}
chroms <- chroms %>% mutate(chrom = trimws(as.character(chrom)))

mises <- readr::read_csv(mis_csv, show_col_types = FALSE, col_names = TRUE)
if (!all(c("chrom", "start", "end", "type") %in% names(mises))) {
  if (ncol(mises) >= 4) {
    names(mises)[1:4] <- c("chrom","start","end","type")
  } else stop("Second CSV must have three columns: chrom,start,end,type (with header).")
}
mises  <- mises  %>% mutate(chrom = trimws(as.character(chrom)))

mises <- mises %>% dplyr::semi_join(chroms, by = "chrom")
mises <- mises %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  mutate(x0 = pmin(start, end), x1 = pmax(start, end)) %>%
  select(chrom, x0, x1, type)
mises <- mises %>% mutate(type = as.factor(type))

chrom_lengths <- chroms %>% select(chrom, length)
mises <- mises %>%
  inner_join(chrom_lengths, by = "chrom") %>%
  mutate(x0 = pmax(0, pmin(x0, length)),
         x1 = pmax(0, pmin(x1, length))) %>%
  filter(x1 > x0)

#min_vis_bp <- max(5e4, 5e-5 * max(chroms$length, na.rm = TRUE))
min_vis_bp <- as.numeric(scale_factor)

mises <- mises %>%
  mutate(len = x1 - x0,
         d   = pmax(len, min_vis_bp),
         mid = (x0 + x1) / 2,
         x0_vis = pmax(0, pmin(length, mid - d/2)),
         x1_vis = pmax(0, pmin(length, mid + d/2)))

chroms$chrom <- factor(chroms$chrom, levels = chroms$chrom)
mises$chrom  <- factor(mises$chrom,  levels = levels(chroms$chrom))

chrom_seg <- chroms %>% transmute(chrom, x0 = 0, x1 = length)
line_mm <- 4

p <- ggplot() +
  geom_segment(data = chrom_seg,
               aes(x = x0, xend = x1, y = chrom, yend = chrom),
               linewidth = line_mm, color = "black", lineend = "round") +
  geom_segment(data = mises,
               aes(x = x0_vis, xend = x1_vis, y = chrom, yend = chrom, colour = type),
               linewidth = line_mm, lineend = "butt") +
  scale_colour_discrete(name = "Misassembly Types")+
  scale_y_discrete(name = NULL, drop = FALSE) +
  scale_x_continuous(name = "Position", expand = expansion(mult = c(0.02, 0.01))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks.y     = element_blank())

max_len <- max(chroms$length, na.rm = TRUE)
n_chr   <- nrow(chroms)
w_in <- max(6, min(16, max_len / 50e6 * 10))
h_in <- max(2.5, min(20, n_chr * 0.35))

ggsave(filename = paste0(out_pref, ".svg"), plot = p, width = w_in, height = h_in, units = "in", device = "svg")
ggsave(filename = paste0(out_pref, ".png"), plot = p, width = w_in, height = h_in, units = "in", dpi = 1200)
