# Green Monster the Runs

library(tidyverse)

zscore_files <- list.files("results", "krsa", full.names = TRUE) |>
  keep(~ str_detect(.x, "all", negate = TRUE)) |>
  set_names(~ str_extract(.x, "results/(run\\d+)", 1))


dataset <- zscore_files |>
  map(~ read_csv(.x, show_col_types = FALSE, col_select = c(Kinase, AvgZ))) |>
  map(unique) |>
  map(~ mutate(.x, rank = row_number(desc(AvgZ)),
               quartile = ntile(desc(AvgZ), 4))) |>
  list_rbind(names_to = "run") |>
  write_csv("results/combined_quartile_ranked_KRSA.csv")

compared <- dataset |>
  select(run, Kinase, quartile) |>
  pivot_wider(names_from = run, values_from = quartile) |>
  mutate(mean_rank = (run01 + run02 + run03) / 3)

top10 <- dataset |>
  filter(Kinase %in% {compared |> slice_min(mean_rank, n = 10) |> pull(Kinase)})

compared |>
  select(-mean_rank) |>
  column_to_rownames("Kinase") |>
  as.matrix() |>
  heatmap()
