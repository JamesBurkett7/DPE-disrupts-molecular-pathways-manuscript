# Create barplot for Kinases of Interest

library(tidyverse)
library(ggprism)

pep_data <- readRDS("datastore/run-all-data_modeled-STK.RDS") |>
  pluck("grouped") |>
  filter(r.seq >= 0.8, str_detect(Peptide, "REF", negate = TRUE), str_detect(Peptide, "^p", negate = TRUE)) |>
  mutate(Group = if_else(Group == "A", "Exposed", "Control"))

mapping <- KRSA::KRSA_Mapping_STK_PamChip_87102_v1

kinases_of_interest <- mapping |> pull(Kinases) |> str_split(" ", simplify = TRUE) |> reduce(c) |> unique() |> keep(~ str_detect(.x, "\\w+")) |> setdiff(c("WNK", "HAL", "STKR", "MSN", "BRSK", "STE11", "NDR"))

combined_data <- mapping |>
  separate_longer_delim(Kinases, delim = " ") |>
  filter(Kinases %in% kinases_of_interest) |>
  rename(Peptide = Substrates, Kinase = Kinases) |>
  inner_join(pep_data) |>
  select(Group, Kinase, Peptide, Slope = slope)

tested <- combined_data |>
  nest(.by = c(Kinase)) |>
  mutate(ttest = map(data, ~ t.test(formula = Slope ~ Group, data = .x)),
         glanced = map(ttest, broom::glance),
         mean_val = map(data, ~ summarise(.x, mean_value = mean(Slope), .by = Group)),
         mean_val = map(mean_val, ~ pivot_wider(.x, names_from = Group, values_from = mean_value, names_prefix = "MeanVal"))) |>
  select(-ttest, -data) |>
  unnest_wider(glanced) |>
  unnest_wider(mean_val) |>
  select(-method, -alternative, -parameter, -estimate1, -estimate2) |>
  mutate(Significant = p.value < 0.05) |>
  pivot_longer(cols = starts_with("MeanVal"), names_to = "Group", values_to = "Mean") |>
  mutate(Group = str_remove(Group, "MeanVal"),
         ErrorMin = Mean - abs(conf.low),
         ErrorMax = Mean + abs(conf.high)) |>
  write_csv("results/t-test_all_kinases.csv")

tested_filtered <- tested |>
  filter(Kinase %in% c("PKCH", "PAKB", "PKN", "ERK")) |>
  write_csv("results/t-test_selected_kinases.csv")

g <- ggplot(tested_filtered, aes(x = Kinase, y = Mean, fill = Group, ymin = ErrorMin, ymax = ErrorMax))

p <- g +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(position = position_dodge(0.8), width = 0.7) +
  scale_fill_prism("viridis") +
  theme_prism()

ggsave("figures/barplot-comparison.svg", plot = p, )
