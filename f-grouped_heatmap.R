# Create a grouped heatmap

library(KRSA)
library(tidyverse)

pep_data <- readRDS("datastore/run-all-data_modeled-STK.RDS") |>
  pluck("grouped") |>
  filter(r.seq >= 0.8, str_detect(Peptide, "REF", negate = TRUE), str_detect(Peptide, "^p", negate = TRUE))

peps <- pep_data |>
  pull(Peptide) |>
  unique()

krsa_heatmap_grouped(pep_data, peps, scale = "row")

dev.off()
pdf(file = "figures/combined_heatmap.pdf")
krsa_heatmap_grouped(pep_data, peps, scale = "row")
dev.off()
