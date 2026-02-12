########################################
### Plot averaged Re and DPP results ###
########################################

library(tidyverse)
library(latex2exp)
library(beastio)
library(readxl)
library(patchwork)

# LEO TODO: Something's up with number of Res in large clusters. Check xml.

# Set theme
my_theme <- theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
theme_set(my_theme)

# Parse logs and remove burnin. ESS > 200 checked with beastiary
dpp <- beastio::readLog(
  filenames = "dpp.log",
  burnin = 0.1,
  as.mcmc = FALSE # return data.frame
)
avg <- beastio::readLog(
  filenames = "re_averaged.log",
  burnin = 0.1,
  as.mcmc = FALSE # return data.frame
)
lge <- c("cluster_168.log", "cluster_161.log", "cluster_37.log") %>%
  lapply(function(x) {
    beastio::readLog(
      filenames = x,
      burnin = 0.1,
      as.mcmc = FALSE # return data.frame
    ) %>%
      mutate(cluster = gsub(x, pattern = "cluster_|\\.log", replacement = ""))
  }) %>%
  bind_rows() %>%
  mutate(cluster = paste("Cluster", cluster))

# Parse metadata and sampling metadata
metadata <- read_excel("sequences/large_clusters.xlsx") %>%
  rename_with(tolower) %>%
  mutate(cluster = gsub(pattern = "_", replacement = " ", cluster)) %>%
  mutate(sampling_time = as.Date(`date of collection`))

# Plots

## Average Re trajetory
oldest <- "2014-01-01" # True value in 2005, but use 2014 to condense plot
youngest <- as.character(max(metadata$sampling_time))

re_intervals <- tibble(
  interval = paste0("ReEpi.", 1:4),
  start = as.Date(c("2022-01-01", "2020-01-01", "2017-01-01", oldest)),
  end = as.Date(c(youngest, "2021-12-31", "2019-12-31", "2016-12-31"))
) %>%
  pivot_longer(
    cols = c(start, end),
    names_to = "date_type",
    values_to = "date"
  )

midpoints <- as.Date(c("2015-07-01", "2018-07-01", "2021-01-01", "2023-01-01"))
labs <- c("Before PrEP", "PrEP rollout", "COVID-19\nlockdown", "Post\nlockdown")

# Colour panels
bg_cols <- alpha(c("white", "grey"), 0.2)
bg_panels <- tibble(
  xmin = re_intervals %>% filter(date_type == "start") %>% pull(date),
  xmax = re_intervals %>% filter(date_type == "end") %>% pull(date),
  ymin = 0,
  ymax = Inf,
  fill = rep(bg_cols, 2)
)

re_traj <- avg %>%
  select(starts_with("ReEpi")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "interval",
    values_to = "Re"
  ) %>%
  group_by(interval) %>%
  summarize(
    mean = mean(Re),
    lower = quantile(Re, 0.025),
    upper = quantile(Re, 0.975),
    .groups = "drop"
  ) %>%
  left_join(re_intervals, by = "interval")

p_re_avg <- ggplot(re_traj, aes(x = date, y = mean)) +
  geom_rect(
    data = bg_panels,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = bg_cols, guide = "none") +
  geom_hline(yintercept = 1, lty = "dashed", color = "grey") +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "purple", alpha = 0.3
  ) +
  ylim(0, 3) +
  geom_line(aes(group = interval), color = "purple", linewidth = 1.2) +
  # Add vertical lines at interval breaks
  geom_vline(
    xintercept = re_intervals %>% filter(date_type == "start") %>% pull(date),
    linetype = "solid", color = "grey"
  ) +
  annotate("text", x = midpoints, y = 2.9, label = labs) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(y = TeX("$R_e$"), x = "Year")
p_re_avg

ggsave(
  plot = p_re_avg, "results/figures/re_avg_trajectory.jpeg",
  width = 6, height = 4, dpi = 600
)

## Large cluster Re trajectories
metadata_large_clusters <- metadata %>%
  filter(cluster %in% c("Cluster 168", "Cluster 161", "Cluster 37"))

oldest <- as.character(min(metadata_large_clusters$sampling_time))
youngest <- as.character(max(metadata_large_clusters$sampling_time))

re_intervals <- tibble(
    interval = c("Re.4", "Re.3", "Re.2", "Re.1"),
    start = as.Date(c("2022-01-01", "2020-01-01", "2017-01-01", oldest)),
    end = as.Date(c(youngest, "2021-12-31", "2019-12-31", "2016-12-31"))
  ) %>%
  pivot_longer(
    cols = c(start, end),
    names_to = "date_type",
    values_to = "date"
  )

# Colour panels
bg_cols <- alpha(c("grey", "white"), 0.2)
bg_panels <- tibble(
  xmin = re_intervals %>% filter(date_type == "start") %>% pull(date),
  xmax = re_intervals %>% filter(date_type == "end") %>% pull(date),
  ymin = 0,
  ymax = Inf,
  fill = rep(bg_cols, 2)
)

# Summarize Re statistics by interval
lge_re_traj <- lge %>%
  pivot_longer(
    cols = starts_with("Re"),
    names_to = "interval",
    values_to = "Re"
  ) %>%
  filter(!(
    interval == "Re.1" & cluster %in% c("Cluster 161", "Cluster 168")
  )) %>%
  group_by(interval, cluster) %>%
  select(starts_with("Re")) %>%
  summarize(
    mean = mean(Re),
    lower = quantile(Re, 0.025, na.rm = TRUE),
    upper = quantile(Re, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(re_intervals, by = "interval", relationship = "many-to-many")

midpoints <- as.Date(c("2017-01-01", "2020-01-01", "2022-01-01"))

p_lge_re <- ggplot() +
  geom_rect(
    data = bg_panels,
    aes(
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax,
      fill = fill
    ),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = bg_cols, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_hline(yintercept = 1, lty = "dashed", color = "grey") +
  geom_vline(
    xintercept = midpoints,
    linetype = "solid", color = "grey40"
  ) +
  ggnewscale::new_scale_fill() +
  geom_ribbon(
    data = lge_re_traj,
    aes(x = date, ymin = lower, ymax = upper, fill = cluster),
    alpha = 0.2
  ) +
  geom_line(
    data = lge_re_traj,
    aes(x = date, y = mean, group = interval, col = cluster),
    linewidth = 1.2
  ) +
  facet_wrap(~cluster, ncol = 1) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  labs(y = TeX("$R_e$"), x = "Year") +
  theme(legend.position = "none")

p_lge_re

ggsave(
  plot = p_lge_re, "results/figures/top3_re.jpeg",
  width = 5, height = 5, dpi = 600
)

## Top 3 sampling prop
p_lge_samp <- lge %>%
  ggplot(aes(x = sampProp.2, fill = cluster)) +
  geom_histogram(bins = 30, alpha = 0.5) +
  geom_vline(
    data = . %>% group_by(cluster) %>% summarize(mean_samp = mean(sampProp.2)),
    aes(xintercept = mean_samp),
    linewidth = 1, linetype = "dashed"
  ) +
  labs(x = "Sampling proportion", y = "Frequency") +
  theme(legend.position = "none") +
  facet_wrap(~cluster, ncol = 1)

p_lge_samp

ggsave(
  plot = p_lge_samp, "results/figures/top_3_p.jpeg",
  width = 5, height = 5, dpi = 600
)

## Top 3 Areas
p_lge_area <- metadata_large_clusters %>%
  group_by(cluster, gws_areas) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(pc = 100 * n / sum(n)) %>%
  mutate(gws_areas = case_when(
    gws_areas == "Sydney (gay postcodes)" ~ "Inner City",
    TRUE ~ gws_areas
  )) %>%
  ggplot(aes(x = gws_areas, y = pc, group = gws_areas, fill = cluster)) +
  geom_bar(
    stat = "identity", position = "dodge",
    alpha = 0.5
  ) +
  geom_text(aes(label = sprintf("%.0f%%", pc)),
    vjust = -0.3,
    size = 3
  ) +
  ylim(0, 77) +
  labs(x = "Area", y = "Percentage of samples (%)") +
  theme(legend.position = "none") +
  facet_wrap(~cluster, ncol = 1)

p_lge_area

ggsave(
  plot = p_lge_area, "results/figures/top3_area_proportion.jpeg",
  width = 5, height = 5, dpi = 600
)

## DPP num Re values
n_re_stats <- dpp %>%
  summarize(
    mean = mean(uniqueReCount),
    lower = quantile(uniqueReCount, 0.025),
    upper = quantile(uniqueReCount, 0.975)
  )

p_dpp_n <- dpp %>%
  select("uniqueReCount") %>%
  ggplot() +
  geom_bar(aes(x = uniqueReCount), fill = "dodgerblue", alpha = 0.5) +
  theme_minimal() +
  scale_x_continuous(breaks = 1:16) +
  labs(x = TeX("Number of $R_e$ values"), y = "Frequency") +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
    label = sprintf(
      "Mean: %.2f\n95%% HPD: [%.0f, %.0f]",
      n_re_stats$mean, n_re_stats$lower, n_re_stats$upper
    ),
    size = 4
  )
ggsave(
  plot = p_dpp_n, "results/figures/dpp_unique_re.jpeg",
  width = 6, height = 4, dpi = 600
)

## Panel plot using patchwork
p_re_avg / (p_lge_re | p_lge_area) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(3.5, 5))

ggsave(
  plot = last_plot(), "results/figures/panel_re.jpeg",
  width = 8, height = 8, dpi = 600
)

ggsave(
  plot = p_lge_re, "results/figures/top3_re.jpeg",
  width = 6, height = 6, dpi = 600
)

## Further plots for presentation
(p_lge_re | p_lge_area) +
  plot_annotation(tag_levels = "A")

ggsave(
  plot = last_plot(), "results/figures/top3_re_area.jpeg",
  width = 8, height = 5, dpi = 600
)

## Cluster timeline with Re Intervals
cluster_timeline <- metadata %>%
  group_by(cluster) %>%
  summarize(
    min_date = min(sampling_time),
    max_date = max(sampling_time),
    .groups = "drop"
  )

# Set cluster factor levels by order of first sample
cluster_order <- cluster_timeline %>%
  arrange(min_date) %>%
  pull(cluster)
metadata$cluster <- factor(metadata$cluster, levels = cluster_order)
cluster_timeline$cluster <- factor(
  cluster_timeline$cluster,
  levels = cluster_order
)

## Timeline plot
oldest <- "2005-01-01" # True value in 2005, but use 2014 to condense plot
youngest <- as.character(max(metadata$sampling_time))

# TODO: Double check interval order
re_intervals <- tibble(
    start = as.Date(c("2022-01-01", "2020-01-01", "2017-01-01", oldest)),
    end = as.Date(c(youngest, "2021-12-31", "2019-12-31", "2016-12-31"))
  ) %>%
  pivot_longer(
    cols = c(start, end),
    names_to = "date_type",
    values_to = "date"
  )

midpoints <- as.Date(c("2010-07-01", "2018-07-01", "2021-01-01", "2023-01-01"))
labs <- c(
  "\nBefore PrEP", "\nPrEP rollout",
  "\nCOVID-19\nlockdown", "\nPost\nlockdown"
)

# Colour panels
bg_cols <- alpha(c("grey", "white"), 0.2)
bg_panels <- tibble(
  xmin = re_intervals %>% filter(date_type == "start") %>% pull(date),
  xmax = re_intervals %>% filter(date_type == "end") %>% pull(date),
  ymin = 0,
  ymax = Inf,
  fill = rep(bg_cols, 2)
)

ggplot() +
  scale_fill_manual(values = bg_cols, guide = "none") +
  geom_segment(
    data = cluster_timeline,
    aes(x = min_date, xend = max_date, y = cluster, yend = cluster),
    color = "grey40", linewidth = 1
  ) +
  geom_point(
    data = metadata,
    aes(x = sampling_time, y = cluster),
    shape = 21, size = 2,
    fill = "dodgerblue", alpha = 0.7
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Sample date", y = "Cluster")

ggsave(
  "results/figures/cluster_timeline_unmarked.jpeg",
  width = 8, height = 5, dpi = 600
)

ggplot() +
  geom_rect(
      data = bg_panels,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      inherit.aes = FALSE
    ) +
  scale_fill_manual(values = bg_cols, guide = "none") +
  geom_vline(
    xintercept = re_intervals %>% filter(date_type == "start") %>% pull(date),
    linetype = "solid", color = "grey"
  ) +
  annotate("text",
    label = labs,
    x = midpoints, y = Inf + 0.5, # Place above the top cluster
    vjust = -0.1, size = 4
  ) +
  geom_segment(
    data = cluster_timeline,
    aes(x = min_date, xend = max_date, y = cluster, yend = cluster),
    color = "grey40", size = 1
  ) +
  geom_point(
    data = metadata,
    aes(x = sampling_time, y = cluster),
    shape = 21, size = 2,
    fill = "dodgerblue", alpha = 0.7
  ) +
  coord_cartesian(clip = 'off') +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Sample date", y = "Cluster") +
  theme(
    plot.margin = margin(50, 5.5, 5.5, 5.5)
  ) # Add right margin for annotation

ggsave(
  "results/figures/cluster_timeline_marked.jpeg",
  width = 9, height = 5.2, dpi = 600
)

## Marked and coloured
ggplot() +
  geom_rect(
      data = bg_panels,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      inherit.aes = FALSE
    ) +
  scale_fill_manual(values = bg_cols, guide = "none") +
  geom_vline(
    xintercept = re_intervals %>% filter(date_type == "start") %>% pull(date),
    linetype = "solid", color = "grey"
  ) +
  annotate("text",
    label = labs,
    x = midpoints, y = Inf + 0.5, # Place above the top cluster
    vjust = -0.1, size = 4
  ) +
  geom_segment(
    data = cluster_timeline,
    aes(x = min_date, xend = max_date, y = cluster, yend = cluster),
    color = "grey40", size = 1
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = metadata,
    aes(x = sampling_time, y = cluster, fill = gws_areas, shape = exposure),
    shape = 21, size = 2,
  alpha = 0.7
  ) +
  coord_cartesian(clip = 'off') +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Sample date", y = "Cluster") +
  theme(
    plot.margin = margin(50, 5.5, 5.5, 5.5)
  ) # Add right margin for annotation

ggsave(
  "results/figures/cluster_timeline_marked_coloured.jpeg",
  width = 9, height = 5.2, dpi = 600
)

## Area over time for samples
cumulative_samples <- metadata %>%
  arrange(sampling_time) %>%
  mutate(gws_areas = case_when(
    gws_areas == "Sydney (gay postcodes)" ~ "Inner City",
    TRUE ~ gws_areas
  )) %>%
  count(gws_areas, sampling_time) %>%
  complete(gws_areas, sampling_time, fill = list(n = 0)) %>%
  arrange(gws_areas, sampling_time) %>%
  group_by(gws_areas) %>%
  mutate(cumulative = cumsum(n))


p_area_flame_stacked <- ggplot() +
  geom_rect(
    data = bg_panels,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = bg_cols, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_area(
    data = cumulative_samples,
    aes(x = sampling_time, y = cumulative / sum(cumulative), fill = gws_areas),
    position = "stack"
  ) +
  scale_fill_viridis_d(name = "Area") +
  labs(x = "Date", y = "Cumulative number of samples") +
  theme(legend.position = "bottom")

ggsave(
  plot = p_area_flame_stacked, "results/figures/area_flame_stacked.jpeg",
  width = 8, height = 5, dpi = 600
)

### Get HPDs for Re ###
re_table <- avg %>%
  select(starts_with("ReEpi")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "interval",
    values_to = "Re"
  ) %>%
  mutate(interval = case_when(
    interval == "ReEpi.1" ~ "Post-lockdown",
    interval == "ReEpi.2" ~ "COVID-19 lockdown",
    interval == "ReEpi.3" ~ "PrEP rollout",
    interval == "ReEpi.4" ~ "Before PrEP"
  )) %>%
  mutate(interval = paste0("Re ", interval)) %>%
  group_by(interval) %>%
  summarize(
    mean = mean(Re),
    lower = quantile(Re, 0.025),
    upper = quantile(Re, 0.975),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .)))

knitr::kable(
  re_table,
  format = "simple",
  col.names = c("Interval", "Mean", "Lower", "Upper")
)
write_tsv(re_table, "results/tables/re_hpd.csv")

### HPS for Re intervals for each large cluster ###
lge_re_table <- lge %>%
  pivot_longer(
    cols = starts_with("Re"),
    names_to = "interval",
    values_to = "Re"
  ) %>%
  filter(!(
    (interval == "Re.1" & cluster %in% c("Cluster 161", "Cluster 168")) |
    (interval %in% c("Re.4", "Re.3") & cluster == "Cluster 37")
  )) %>%
  # Note that interval ordering is reverse compared to in avg_re
  mutate(interval = case_when(
    interval == "Re.4" ~ "Post-lockdown",
    interval == "Re.3" ~ "COVID-19 lockdown",
    interval == "Re.2" ~ "PrEP rollout",
    interval == "Re.1" ~ "Before PrEP"
  )) %>%
  group_by(cluster, interval) %>%
  summarize(
    mean = mean(Re, na.rm = TRUE),
    lower = quantile(Re, 0.025, na.rm = TRUE),
    upper = quantile(Re, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .)))
write_tsv(lge_re_table, "results/tables/lge_clusters_re_hpd.tsv")

### Get HPDs for DPP number of Re values
dpp_table <- dpp %>%
  summarize(
    mean = mean(uniqueReCount),
    lower = quantile(uniqueReCount, 0.025),
    upper = quantile(uniqueReCount, 0.975)
  ) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .)))
write_tsv(dpp_table, "results/tables/dpp_n_re_hpd.tsv")

## DPP Heatmap ##
dpp_long <- dpp %>%
  select(starts_with("Re_Cluster_")) %>%
  mutate(state = row_number()) %>%
  pivot_longer(cols = -state, names_to = "cluster", values_to = "Re")

similarity <- dpp_long %>%
  inner_join(dpp_long, by = c("state", "Re")) %>%
  filter(cluster.x != cluster.y) %>%
  count(cluster.x, cluster.y, name = "matches") %>%
  mutate(
    cluster.x = gsub("Re_Cluster_", "", cluster.x),
    cluster.y = gsub("Re_Cluster_", "", cluster.y),
    prop_matches = matches / max(dpp_long$state)
  )
  # Calculate cluster order based on total prop_matches
  cluster_order <- similarity %>%
    group_by(cluster.x) %>%
    summarize(total_matches = sum(prop_matches)) %>%
    arrange(desc(total_matches)) %>%
    pull(cluster.x)

  # Convert to factor with sorted levels
  similarity <- similarity %>%
    mutate(
      cluster.x = factor(cluster.x, levels = cluster_order),
      cluster.y = factor(cluster.y, levels = cluster_order)
    )

  ggplot(similarity, aes(x = cluster.x, y = cluster.y, fill = prop_matches)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Number of\nmatching states") +
    labs(x = "Cluster", y = "Cluster") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      text = element_text(size = 12)
    ) +
    coord_fixed()

## DPP Heatmap with marginal stacked barplot ##
# Remove state where all are equal.
dpp_long <- dpp %>%
  filter(uniqueReCount > 1) %>%
  select(starts_with("Re_Cluster_")) %>%
  mutate(state = row_number()) %>%
  pivot_longer(cols = -state, names_to = "cluster", values_to = "Re")

similarity <- dpp_long %>%
  inner_join(dpp_long, by = c("state", "Re")) %>%
  filter(cluster.x != cluster.y) %>%
  count(cluster.x, cluster.y, name = "matches") %>%
  mutate(
    cluster.x = gsub("Re_Cluster_", "", cluster.x),
    cluster.y = gsub("Re_Cluster_", "", cluster.y),
    prop_matches = matches / max(dpp_long$state)
  )
  # Calculate cluster order based on total prop_matches
  cluster_order <- similarity %>%
    group_by(cluster.x) %>%
    summarize(total_matches = sum(prop_matches)) %>%
    arrange(desc(total_matches)) %>%
    pull(cluster.x)

  # Convert to factor with sorted levels
  similarity <- similarity %>%
    mutate(
      cluster.x = factor(cluster.x, levels = cluster_order),
      cluster.y = factor(cluster.y, levels = cluster_order)
    )


# Calculate area proportions per cluster
area_prop <- metadata %>%
  mutate(gws_areas = case_when(
    gws_areas == "Sydney (gay postcodes)" ~ "Inner City",
    TRUE ~ gws_areas
  )) %>%
  mutate(cluster_id = gsub("Cluster ", "", cluster)) %>%
  group_by(cluster_id, gws_areas) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(cluster_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    cluster_id = factor(cluster_id, levels = cluster_order),
    gws_areas = factor(gws_areas, levels = c("GWS", "Inner City", "Other Sydney", "Rest of NSW"))
  )

# Main heatmap
p_heatmap <- ggplot(similarity, aes(x = cluster.x, y = cluster.y, fill = prop_matches)) +
  geom_tile() +
  scale_fill_viridis_c(name = TeX("Proportion matching $R_0$")) +
  labs(x = NULL, y = "Cluster") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.position = "right"
  ) +
  theme(
    text = element_text(size = 14),
    legend.position = "top"
  )

# Marginal stacked barplot for area proportions
area_cols <- c(
  "GWS" = "red",
  "Inner City" = "dodgerblue",
  "Other Sydney" = "#66C2A5",
  "Rest of NSW" = "#FEE08B"
)

p_marginal <- ggplot(area_prop, aes(x = cluster_id, y = prop, fill = gws_areas)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = area_cols, name = "Area") +
  labs(x = "Cluster", y = "Proportion") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  theme(text = element_text(size = 14))

# Combine with patchwork
p_heatmap_marginal <- p_heatmap / p_marginal +
  plot_layout(heights = c(5, 1), axes = "collect") &
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

p_heatmap_marginal

ggsave(
  plot = p_heatmap_marginal,
  "results/figures/dpp_heatmap_area.jpeg",
  dpi = 600
)
