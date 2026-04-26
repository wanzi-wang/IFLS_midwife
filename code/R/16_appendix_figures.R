# =====================================================================
# 16_appendix_figures.R — Appendix-only figures.
#
# Outputs:
#   paper/figures/midwife_rollout_map.pdf — three-panel province-level
#     choropleth showing the share of IFLS communities with a village
#     midwife present in 1990, 1994, and 1999.
#
# IFLS-13 province codes (BPS) → rnaturalearth name mapping is built
# below. Non-IFLS provinces are drawn in light grey for geographic
# context.
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthhires)
  library(patchwork)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))

log_stage("stage16", "begin")

frame       <- read_intermediate("stage06", "analysis_sample")
figures_dir <- here::here("paper", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1. Province-by-snapshot-year share of IFLS communities with midwife.
# ---------------------------------------------------------------------
prov_bps <- tribble(
  ~bps_code, ~ne_name,
  12L, "Sumatera Utara",
  13L, "Sumatera Barat",
  16L, "Sumatera Selatan",
  18L, "Lampung",
  31L, "Jakarta Raya",
  32L, "Jawa Barat",
  33L, "Jawa Tengah",
  34L, "Yogyakarta",
  35L, "Jawa Timur",
  51L, "Bali",
  52L, "Nusa Tenggara Barat",
  63L, "Kalimantan Selatan",
  73L, "Sulawesi Selatan"
)

comm_level <- frame |>
  filter(primary_sample == 1) |>
  distinct(commid_birth_legacy, province, start_sar_legacy) |>
  filter(province %in% prov_bps$bps_code)

snapshot_years <- c(1990, 1994, 1999)

shares <- bind_rows(lapply(snapshot_years, function(y) {
  comm_level |>
    group_by(province) |>
    summarise(
      n_comm  = dplyr::n(),
      n_treat = sum(!is.na(start_sar_legacy) & start_sar_legacy <= y),
      share   = 100 * n_treat / n_comm,
      .groups = "drop"
    ) |>
    mutate(year = y)
})) |>
  left_join(prov_bps, by = c("province" = "bps_code"))

# ---------------------------------------------------------------------
# 2. Province polygons.
# ---------------------------------------------------------------------
prov_sf <- ne_states(country = "Indonesia", returnclass = "sf")

prov_sf_with_share <- bind_rows(lapply(snapshot_years, function(y) {
  prov_sf |>
    left_join(shares |> filter(year == y),
              by = c("name" = "ne_name")) |>
    mutate(year = y)
}))

# ---------------------------------------------------------------------
# 3. Plot — one panel per snapshot year, shared color scale.
# ---------------------------------------------------------------------
make_panel <- function(yy) {
  d <- prov_sf_with_share |> filter(year == yy)
  ggplot(d) +
    geom_sf(aes(fill = share), color = "black", linewidth = 0.15) +
    scale_fill_gradient(
      name   = "% IFLS communities",
      low    = "#D6F0E5", high = "#08306B",
      na.value = "white", limits = c(0, 100)
    ) +
    labs(title = sprintf("Share of IFLS communities with midwife in %d", yy)) +
    coord_sf(xlim = c(94, 142), ylim = c(-11, 6), expand = FALSE) +
    theme_classic(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.position  = "bottom",
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(1.2, "cm"),
      plot.title       = element_text(size = 10)
    )
}

panels <- lapply(snapshot_years, make_panel)
p_combined <- panels[[1]] / panels[[2]] / panels[[3]] +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "midwife_rollout_map.pdf"),
       plot = p_combined, width = 7.5, height = 9.0, units = "in",
       device = cairo_pdf, bg = "white")

log_stage("stage16", sprintf(
  "Midwife rollout map written for snapshot years %s.",
  paste(snapshot_years, collapse = ", ")))
log_stage("stage16", "done")
