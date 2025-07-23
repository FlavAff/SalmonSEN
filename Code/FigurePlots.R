setwd("~/Documents/McGill/PhD/Chapter 3/Code/")

library(tidyverse)
library(posterior)
library(patchwork)
library(ggh4x)


Chindat <- SES |> filter(species == "Chinook")
Cohodat <- SES |> filter(species == "Coho")


#Get the right model and data
dat <- SES |> filter(species == "Coho")
Mod <- CohoPreds

predictions_long <- Mod |>
  gather_rvars(
    log_Erec_rep[region.n, year.n], log_Crec[region.n, year.n],
    log_E_rep[region.n, year.n], log_E[region.n, year.n],
    logit_hrec[region.n, year.n],
    logit_h[region.n, year.n]
  )

predictions_wide <- pivot_wider(predictions_long, names_from = .variable, values_from = .value)
plot_data <- predictions_wide |>
  left_join(dat, by = c("region.n","year.n"))

CohoPlot <- plot_data
ChinookPlot <- plot_data

plot_data <- rbind(ChinookPlot, CohoPlot)

# effort_dat <- plot_data |> 
#   group_by(year, species) |> 
#   summarise(comm_rvar = rvar_sum(log_E_rep, na.rm = TRUE),
#             rec_rvar = rvar_sum(log_Erec_rep, na.rm = TRUE),
#             comm_effort = sum(total_effort, na.rm = T),
#             rec_effort = sum(total_rec.effort/100, na.rm = T),
#             .groups = 'drop')
# 
# effort_dat |>
#   ggplot(aes(x = year, dist = comm_rvar)) +
#   stat_lineribbon() +
#   facet_wrap(. ~ species) +
#   geom_point(data = effort_dat, aes(x = year, y = log(comm_effort)),
#              colour = "blue", alpha = 0.5, inherit.aes = FALSE) +
#   scale_fill_brewer(palette = "Oranges", direction = -1) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "Year", y = "Predicted log abundance")
# 
# effort_dat |>
#   ggplot(aes(x = year, dist = rec_rvar)) +
#   stat_lineribbon() +
#   facet_wrap(. ~ species) +
#   geom_point(data = Cohodat, aes(x = year, y = log(total_rec.effort/100)),
#              colour = "pink", alpha = 0.5, inherit.aes = FALSE) +
#   scale_fill_brewer(palette = "BuGn", direction = -1) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "Year", y = "Predicted log abundance")


#Average effort for both species across regions
effort_dat2 <- plot_data |> 
  group_by(year, region) |> 
  summarise(comm_rvar = rvar_mean(log_E_rep, na.rm = TRUE),
            rec_rvar = rvar_mean(log_Erec_rep, na.rm = TRUE),
            .groups = 'drop')

Ec_p <- effort_dat2 |>
  ggplot(aes(x = year, dist = comm_rvar)) +
  stat_lineribbon() +
  facet_wrap(. ~ region) +
  scale_fill_brewer(palette = "Oranges", direction = -1) +
  geom_point(data = Chindat, aes(x = year, y = log(total_effort)),
             colour = "blue", alpha = 0.5, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Year", y = "Average commercial effort (log)")+
  theme_minimal(base_size = 18)
ggsave(filename = "../Results/HPCfits/AvgEffortC.png", plot = Ec_p, width = 9, height = 6, dpi = 300)

Er_p <- effort_dat2 |>
  ggplot(aes(x = year, dist = rec_rvar)) +
  stat_lineribbon() +
  facet_wrap(. ~ region) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  geom_point(data = Chindat, aes(x = year, y = log(total_rec.effort/100)),
             colour = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Year", y = "Average recreational effort (log)")+
  theme_minimal(base_size = 18)
ggsave(filename = "../Results/HPCfits/AvgEffortR.png", plot = Er_p, width = 9, height = 6, dpi = 300)



# Summarize harvest rates for each year and species
harvest_rate_summary <- plot_data %>%
  group_by(species, year) %>%
  summarise(
    # Get the mean harvest rate across locations for each draw
    logit_h_comm = rvar_mean(logit_h),
    logit_h_rec = rvar_mean(logit_hrec),
    .groups = "drop"
  )

harvest_rate_summary$h_comm <- exp(harvest_rate_summary$logit_h_comm)/(1+exp(harvest_rate_summary$logit_h_comm)+exp(harvest_rate_summary$logit_h_rec))
harvest_rate_summary$h_rec <- exp(harvest_rate_summary$logit_h_rec)/(1+exp(harvest_rate_summary$logit_h_comm)+exp(harvest_rate_summary$logit_h_rec))

# Sample draws from an rvar
sample_draws <- function(rv, n = 200) {
  draws <- draws_of(rv)
  if (length(draws) <= n) return(draws)
  sample(draws, n)
}

# Set number of draws you want per rvar (per row)
n_draws <- 1000

# Apply to your dataframe
harvest_rate <- harvest_rate_summary %>%
  transmute(
    id = row_number(),
    species = species,                                
    h_draws = map(h_comm, sample_draws, n = n_draws),
    hrec_draws = map(h_rec, sample_draws, n = n_draws)
  ) %>%
  unnest(c(h_draws, hrec_draws), names_repair = "minimal")


h_plot <- ggplot(harvest_rate, aes(x = h_draws, y = hrec_draws)) +
  # Points: light colors
  geom_point(aes(color = species), alpha = 0.3) +
  
  # Lines: manually set dark colors by species
  geom_smooth(
    aes(group = species),  # ensure one line per species
    method = "lm", se = TRUE, size = 1.2,
    color = NA  # remove default mapping
  ) +
  # Add each species line separately with custom color
  geom_smooth(
    data = subset(harvest_rate, species == "Coho"),
    aes(x = h_draws, y = hrec_draws),
    method = "lm", se = TRUE, color = "darkblue", size = 1.2
  ) +
  geom_smooth(
    data = subset(harvest_rate, species == "Chinook"),
    aes(x = h_draws, y = hrec_draws),
    method = "lm", se = TRUE, color = "darkred", size = 1.2
  ) +
  
  # Points color scale: light versions
  scale_color_manual(values = c(
    "Coho" = "#87CEFA",    # light blue
    "Chinook" = "#F08080"  # light red
  )) +
  
  theme_minimal() +
  labs(
    x = "Commercial harvest rate",
    y = "Recreational harvest rate",
    color = "Species"
  ) +
  theme_minimal(base_size = 18)

ggsave(filename = "../Results/HPCfits/hcomparison.png", plot = h_plot, width = 6, height = 8, dpi = 300)



# Get the catchability plots done

# Function to extract, exponentiate, and prepare parameters from a fit object
process_params <- function(fit, species_id) {
  fit$draws(variables = c("log_q", "log_qrec"), format = "df") %>%
    # Exponentiate all parameter columns
    mutate(across(everything(), exp)) %>%
    # Rename columns to remove "log_" prefix
    rename_with(~str_remove(., "log_"), starts_with("log_")) %>%
    # Pivot to a long format suitable for plotting
    pivot_longer(
      cols = everything(),
      names_to = c("parameter", "region"),
      names_pattern = "([a-zA-Z_]+)\\[(\\d+)\\]",
      values_to = "value"
    ) %>%
    # Add a column identifying the species
    mutate(species = species_id)
}

# Process the data for each species
params_species1 <- process_params(ChinookPreds, "Chinook")
params_species2 <- process_params(CohoPreds, "Coho")

# Combine into a single dataframe
all_params <- bind_rows(params_species1, params_species2)

# Calculate median and 95% credible intervals for plotting
plot_summary <- all_params %>%
  group_by(species, parameter, region) %>%
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )
plot_summary <- plot_summary |> filter(region != "NA")

# Define the new labels as a named vector
# The names ("1", "2", etc.) must match the existing values in the 'region' column.
region_labels <- c(
  "1" = "CC",
  "2" = "Fraser",
  "3" = "HG",
  "4" = "Nass",
  "5" = "Skeena",
  "6" = "VIMI"
)

# Generate the plot with renamed x-axis labels
q_plot <- ggplot(plot_summary, aes(x = region, y = median, color = parameter)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.5
  ) +
  geom_point(
    position = position_dodge(width = 0.7),
    size = 3
  ) +
  facet_wrap(~species, scales = "free_y") +
  ggh4x::facetted_pos_scales(
    y = list(
      species == "Coho" ~ scale_y_continuous(limits = c(NA, 0.016))
    )
  ) +
  
  # --- Add this line to rename the x-axis ticks ---
  scale_x_discrete(labels = region_labels) +
  # ------------------------------------------------

labs(
  x = "Region",
  y = "Catchability",
  color = "Parameter"
) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  scale_color_brewer(palette = "Dark2",
                     labels = c("q" = expression(q[ ]), "qrec" = expression(q[r])))

ggsave(filename = "../Results/HPCfits/qcomparison.png", plot = q_plot, width = 8, height = 6, dpi = 300)
