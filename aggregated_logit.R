library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(fixest)
library(broom)
library(tibble)
library(stringr)
library(tidyr)
library(ggplot2)


setwd("C:/Users/LENOVO/Documents/GitHub/mixedlogit_uk")

df <- read_excel("df_final.xlsx")

# Identify the 77 subgroup count columns
groups <- names(df) %>%
  str_subset(
    "^(white|mixed|asian|black)_.*_sc[1-7]$"
  )

# Convert wide counts to long format
df_choice <- df %>%
  pivot_longer(
    cols = all_of(groups),
    names_to = "group",
    values_to = "n"
  ) %>%
  mutate(
    ethnic_group = str_remove(group, "_sc[1-7]$"),
    ethnic_group = str_remove(
      ethnic_group,
      "^(white|mixed|asian|black)_"
    ),
    sc_group = str_extract(group, "sc[1-7]"),
    year = as.integer(year),
    group_year = paste(group, year, sep = "_")
  )

# Retrieve the matching concentration variables
ethnic_cols <- paste0(
  "frac_total_",
  df_choice$ethnic_group
)

class_cols <- paste0(
  "frac_total_",
  df_choice$sc_group
)

df_choice$frac_homoe <-
  as.numeric(
    as.matrix(df_choice)[
      cbind(
        seq_len(nrow(df_choice)),
        match(ethnic_cols, names(df_choice))
      )
    ]
  )

df_choice$frac_homos <-
  as.numeric(
    as.matrix(df_choice)[
      cbind(
        seq_len(nrow(df_choice)),
        match(class_cols, names(df_choice))
      )
    ]
  )

df_choice %>%
  group_by(group_year) %>%
  summarise(
    n_boroughs = n_distinct(borough),
    total_residents = sum(n),
    .groups = "drop"
  ) %>%
  print(n = 20)

stopifnot(
  !anyNA(df_choice$frac_homoe),
  !anyNA(df_choice$frac_homos),
  !anyNA(df_choice$n),
  all(df_choice$pop_borough > 0),
  all(df_choice$n >= 0),
  all(df_choice$frac_homoe >= 0 &
        df_choice$frac_homoe <= 1),
  all(df_choice$frac_homos >= 0 &
        df_choice$frac_homos <= 1)
)

mod_A <- fepois(
  n ~
    i(group_year, frac_homoe) +
    i(group_year, frac_homos) |
    group_year,
  data = df_choice,
  fixef.rm = "none"
)

# Model B: population-adjusted model
mod_B <- fepois(
  n ~
    i(group_year, frac_homoe) +
    i(group_year, frac_homos) |
    group_year,
  offset = ~ log(pop_borough),
  data = df_choice,
  fixef.rm = "none"
)

summary(mod_B)

results <- tidy(mod_B) %>%
  filter(
    grepl("frac_homoe|frac_homos", term)
  ) %>%
  mutate(
    dimension = if_else(
      grepl("frac_homoe", term),
      "Ethnic",
      "Social class"
    ),
    OR_10pp = exp(estimate * 0.10)
  )

head(results)

identification_check <- df_choice %>%
  group_by(group, year) %>%
  summarise(
    total_n = sum(n),
    n_boroughs = n_distinct(borough),
    sd_ethnic = sd(frac_homoe),
    sd_class = sd(frac_homos),
    correlation = cor(
      frac_homoe,
      frac_homos
    ),
    .groups = "drop"
  )

summary(identification_check)

results <- tibble(
  term = names(coef(mod_B)),
  estimate = unname(coef(mod_B))
) %>%
  filter(str_detect(term, "frac_homoe|frac_homos")) %>%
  mutate(
    dimension = if_else(
      str_detect(term, "frac_homoe$"),
      "Ethnic",
      "Social class"
    ),
    group_year = str_remove(
      term,
      "^group_year::"
    ) %>%
      str_remove(":frac_homoe$") %>%
      str_remove(":frac_homos$"),
    year = as.integer(str_extract(group_year, "\\d{4}$")),
    group = str_remove(group_year, "_\\d{4}$")
  ) %>%
  select(group, year, dimension, estimate)

# Calculate IQRs within each subgroup and wave
iqr_data <- df_choice %>%
  group_by(group, year) %>%
  summarise(
    Ethnic = IQR(frac_homoe, na.rm = TRUE),
    `Social class` = IQR(frac_homos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Ethnic, `Social class`),
    names_to = "dimension",
    values_to = "IQR"
  )

# Join estimated coefficients to the IQRs
results_iqr <- results %>%
  mutate(year = as.integer(year)) %>%
  left_join(
    iqr_data,
    by = c("group", "year", "dimension")
  ) %>%
  mutate(
    effect_IQR = estimate * IQR,
    OR_IQR = exp(effect_IQR)
  )

results_iqr %>%
  filter(group %in% c(
    "asian_indian_sc1",
    "asian_chinese_sc1",
    "asian_indian_sc7"
  )) %>%
  arrange(group, year, dimension)

ggplot(
  results_iqr,
  aes(
    x = factor(year),
    y = effect_IQR,
    color = dimension,
    group = dimension
  )
) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~group, scales = "free_y") +
  labs(
    x = "Census wave",
    y = "IQR-adjusted log relative odds",
    color = "Dimension"
  ) +
  theme_bw()

ggplot(
  results_iqr,
  aes(
    x = factor(year),
    y = group,
    fill = effect_IQR
  )
) +
  geom_tile() +
  facet_wrap(~dimension, nrow = 1) +
  scale_fill_gradient2(midpoint = 0) +
  labs(
    x = "Census wave",
    y = "Ethnic–social-class subgroup",
    fill = "IQR effect"
  ) +
  theme_minimal()

comparison <- results_iqr %>%
  select(group, year, dimension, effect_IQR) %>%
  pivot_wider(
    names_from = dimension,
    values_from = effect_IQR
  ) %>%
  mutate(
    difference = Ethnic - `Social class`
  )

ggplot(
  comparison,
  aes(
    x = factor(year),
    y = reorder(group, difference),
    fill = difference
  )
) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0) +
  labs(
    x = "Census wave",
    y = "Subgroup",
    fill = "Ethnic − class"
  ) +
  theme_minimal()

df_prob <- df_choice %>%
  mutate(
    predicted_n = as.numeric(
      predict(mod_B, newdata = df_choice, type = "response")
    )
  ) %>%
  group_by(group, year) %>%
  mutate(
    predicted_prob = predicted_n / sum(predicted_n),
    observed_prob = n / sum(n)
  ) %>%
  ungroup()

model_fit <- df_prob %>%
  group_by(group, year) %>%
  summarise(
    MAE = mean(
      abs(observed_prob - predicted_prob)
    ),
    
    RMSE = sqrt(mean(
      (observed_prob - predicted_prob)^2
    )),
    
    TVD = 0.5 * sum(
      abs(observed_prob - predicted_prob)
    ),
    
    .groups = "drop"
  )

##

coef_df <- tibble(
  term = names(coef(mod_B)),
  estimate = unname(coef(mod_B))
) %>%
  filter(str_detect(term, "frac_homoe|frac_homos")) %>%
  mutate(
    dimension = if_else(
      str_detect(term, ":frac_homoe$"),
      "beta_E",
      "beta_S"
    ),
    group_year = term %>%
      str_remove("^group_year::") %>%
      str_remove(":frac_homoe$") %>%
      str_remove(":frac_homos$")
  ) %>%
  select(group_year, dimension, estimate) %>%
  pivot_wider(
    names_from = dimension,
    values_from = estimate
  )

g <- "white_ewsnib_sc1"
t <- 2001

b <- coef_df %>%
  filter(group_year == paste0(g, "_", t))

beta_E <- b$beta_E
beta_S <- b$beta_S

print(b)

d <- df_choice %>%
  filter(group == g, year == t) %>%
  arrange(borough)

V <- with(
  d,
  log(pop_borough) +
    beta_E * frac_homoe +
    beta_S * frac_homos
)


prob_curve <- function(variable, values) {
  
  bind_rows(lapply(seq_len(nrow(d)), function(j) {
    
    bind_rows(lapply(values, function(x) {
      
      V_new <- V
      
      if (variable == "Ethnic") {
        
        V_new[j] <- log(d$pop_borough[j]) +
          beta_E * x +
          beta_S * d$frac_homos[j]
        
      } else if (variable == "Social class") {
        
        V_new[j] <- log(d$pop_borough[j]) +
          beta_E * d$frac_homoe[j] +
          beta_S * x
        
      } else {
        stop("Use 'Ethnic' or 'Social class'")
      }
      
      # Multinomial probabilities across all 33 boroughs
      w <- exp(V_new - max(V_new))
      p <- w / sum(w)
      
      data.frame(
        borough = d$borough[j],
        concentration = x,
        probability = p[j],
        dimension = variable
      )
    }))
  }))
}

curve_E <- prob_curve(
  variable = "Ethnic",
  values = seq(
    min(d$frac_homoe),
    max(d$frac_homoe),
    length.out = 100
  )
)

curve_S <- prob_curve(
  variable = "Social class",
  values = seq(
    min(d$frac_homos),
    max(d$frac_homos),
    length.out = 100
  )
)

curves_avg <- bind_rows(curve_E, curve_S) %>%
  group_by(dimension, concentration) %>%
  summarise(
    mean_probability = mean(probability),
    .groups = "drop"
  )

ggplot(
  curves_avg,
  aes(
    x = concentration,
    y = mean_probability,
    color = dimension
  )
) +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Indian SC1 — 2021",
    subtitle = "33-alternative conditional logit",
    x = "Neighbourhood concentration",
    y = "Mean predicted selection probability",
    color = "Dimension"
  ) +
  theme_bw()

# margin

df_ame <- df_prob %>%
  left_join(
    coef_df,
    by = "group_year"
  ) %>%
  mutate(
    # Marginal effects of ethnic concentration
    ME_ethnic = beta_E *
      predicted_prob * (1 - predicted_prob),
    
    # Marginal effects of social-class concentration
    ME_class = beta_S *
      predicted_prob * (1 - predicted_prob)
  ) %>%
  group_by(group, year) %>%
  summarise(
    AME_ethnic = mean(ME_ethnic),
    AME_class = mean(ME_class),
    
    # Effects of a 1 percentage-point increase
    AME_ethnic_1pp = AME_ethnic * 0.01,
    AME_class_1pp = AME_class * 0.01,
    
    .groups = "drop"
  )

df_ame %>%
  filter(
    group == "asian_indian_sc1",
    year == 2021
  )

df_ame %>%
  filter(group == "asian_indian_sc1") %>%
  select(
    year,
    AME_ethnic_1pp,
    AME_class_1pp
  ) %>%
  pivot_longer(
    cols = starts_with("AME_"),
    names_to = "dimension",
    values_to = "AME"
  ) %>%
  ggplot(aes(
    x = factor(year),
    y = AME,
    fill = dimension
  )) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Indian SC1 — Average marginal effects",
    x = "Year",
    y = "Change in selection probability (percentage points)",
    fill = "Dimension"
  ) +
  theme_bw()
##
# Select subgroup and year
g <- "asian_chinese_sc1"
year_selected <- 2021

# Extract estimated coefficients
b <- coef_df %>%
  filter(group_year == paste0(g, "_", year_selected))

beta_E <- b$beta_E
beta_S <- b$beta_S

# Original choice set
d <- df_choice %>%
  filter(
    group == g,
    year == year_selected
  ) %>%
  arrange(borough)

# Baseline utilities
V <- with(
  d,
  log(pop_borough) +
    beta_E * frac_homoe +
    beta_S * frac_homos
)

# Simulate probabilities
simulate_prob <- function(variable) {
  
  concentrations <- seq(0, 1, by = 0.01)
  
  bind_rows(lapply(seq_len(nrow(d)), function(j) {
    
    bind_rows(lapply(concentrations, function(x) {
      
      V_new <- V
      
      if (variable == "Ethnic") {
        
        V_new[j] <- log(d$pop_borough[j]) +
          beta_E * x +
          beta_S * d$frac_homos[j]
        
      } else if (variable == "Social class") {
        
        V_new[j] <- log(d$pop_borough[j]) +
          beta_E * d$frac_homoe[j] +
          beta_S * x
      }
      
      # Multinomial probabilities
      w <- exp(V_new - max(V_new))
      P <- w / sum(w)
      
      tibble(
        borough = d$borough[j],
        concentration = x,
        probability = P[j],
        dimension = variable
      )
    }))
  }))
}

# Run simulations
simulated_prob <- bind_rows(
  simulate_prob("Ethnic"),
  simulate_prob("Social class")
)

# Average over the 33 focal borough scenarios
simulated_mean <- simulated_prob %>%
  group_by(dimension, concentration) %>%
  summarise(
    probability = mean(probability),
    .groups = "drop"
  )

ggplot(
  simulated_mean,
  aes(
    x = concentration,
    y = probability,
    color = dimension
  )
) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.1),
    labels = scales::percent
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Indian SC1 — 2021",
    subtitle = "Simulated multinomial residential selection",
    x = "Neighbourhood concentration",
    y = "Mean predicted selection probability",
    color = "Dimension"
  ) +
  theme_bw()


