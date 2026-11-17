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
library(purrr)
library(writexl)

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

# the fraction is put on the scale 0-100 to facilitate reading of log-odds and odds ratio
df_choice <- df_choice %>%
  mutate(
    frac_homoe = frac_homoe * 100,
    frac_homos = frac_homos * 100
  )

df_choice %>%
  group_by(group_year) %>%
  summarise(
    n_boroughs = n_distinct(borough),
    total_residents = sum(n),
    .groups = "drop"
  ) %>%
  print(n = 20)

# sensitivity test (ignore)
stopifnot(
  !anyNA(df_choice$frac_homoe),
  !anyNA(df_choice$frac_homos),
  !anyNA(df_choice$n),
  all(df_choice$pop_borough > 0),
  all(df_choice$n >= 0),
  all(df_choice$frac_homoe >= 0 &
        df_choice$frac_homoe <= 100),
  all(df_choice$frac_homos >= 0 &
        df_choice$frac_homos <= 100)
)

# aggregated conditional model

# only fractions ses and ethnicity
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

# for the quadratic model
df_choice <- df_choice %>%
  mutate(
    frac_homoe_sq = frac_homoe^2,
    frac_homos_sq = frac_homos^2
  )

mod_B_quad <- fepois(
  n ~
    i(group_year, frac_homoe) +
    i(group_year, frac_homoe_sq) +
    i(group_year, frac_homos) +
    i(group_year, frac_homos_sq) |
    group_year,
  offset = ~ log(pop_borough),
  data = df_choice,
  fixef.rm = "none"
)


# comparison of models

dfcomp1 <- mod_B
dfcomp2 <- mod_B_quad
namedfcomp1 <- "mod_B"
namedfcomp2 <- "mod_B_quad"

model_comparison <- tibble(
 Model = c(namedfcomp1,namedfcomp2),
  LogLik = c(
    as.numeric(logLik(dfcomp1)),
    as.numeric(logLik(dfcomp2))
  ),
  AIC = c(AIC(dfcomp1), AIC(dfcomp2)),
  BIC = c(BIC(dfcomp1), BIC(dfcomp2))
)

model_comparison
write_xlsx(model_comparison, paste0("model_comparison/model_comparison","_",namedfcomp1,"_",namedfcomp2,".xlsx"))

# predicted probabilities by the model pred_
# observed probabilities from count number over the total count by groupyear
# predict() gives the number of residential counts, then normalized by population
# to get probabilities

df_compare <- df_choice %>%
  mutate(
    pred_1 = as.numeric(
      predict(dfcomp1, type = "response")
    ),
    pred_2 = as.numeric(
      predict(dfcomp2, type = "response")
    )
  ) %>%
  group_by(group, year) %>%
  mutate(
    observed_prob = n / sum(n),
    
    prob_1 = pred_1 / sum(pred_1),
    prob_2 = pred_2 / sum(pred_2)
  ) %>%
  ungroup()

model_fit <- df_compare %>%
  group_by(group, year) %>%
  summarise(
    TVD_1 = 0.5 * sum(abs(observed_prob - prob_1)),
    TVD_2 = 0.5 * sum(abs(observed_prob - prob_2)),
    
    RMSE_1 = sqrt(mean((observed_prob - prob_1)^2)),
    RMSE_2 = sqrt(mean((observed_prob - prob_2)^2)),
    
    MAE_1 = mean(abs(observed_prob - prob_1)),
    MAE_2 = mean(abs(observed_prob - prob_2)),
    
    .groups = "drop"
  )

model_fit %>%
  summarise(
    mean_TVD_1 = mean(TVD_1),
    mean_TVD_2 = mean(TVD_2),
    
    mean_RMSE_1 = mean(RMSE_1),
    mean_RMSE_2 = mean(RMSE_2),
    
    better_2_TVD = sum(TVD_2 < TVD_1),
    better_1_TVD = sum(TVD_1 < TVD_2)
  )  %>%
  write_xlsx(paste0("model_comparison/model_fit","_",namedfcomp1,"_",namedfcomp2,".xlsx"))

#

mod_B_quad %>%
  coeftable() %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  mutate(
    significance = case_when(
      `Pr(>|z|)` < 0.001 ~ "***",
      `Pr(>|z|)` < 0.01  ~ "**",
      `Pr(>|z|)` < 0.05  ~ "*",
      `Pr(>|z|)` < 0.1   ~ ".",
      TRUE ~ ""
    )
  ) %>%
  write_xlsx("summary_modelsquad.xlsx")
# *** p < 0.001; ** p < 0.01; * p < 0.05; . p < 0.1

# for linear coefficients SC1 - SC7 comparison

coef_df <- tibble(
  term = names(coef(mod_B_quad)),
  estimate = unname(coef(mod_B_quad))
) %>%
  filter(str_detect(term, ":(frac_homoe|frac_homos)$")) %>%
  mutate(
    dimension = if_else(
      str_detect(term, ":frac_homoe$"),
      "beta_E",
      "beta_S"
    ),
    group_year = term %>%
      str_remove("^group_year::") %>%
      str_remove(":(frac_homoe|frac_homos)$")
  ) %>%
  select(group_year, dimension, estimate) %>%
  pivot_wider(
    names_from = dimension,
    values_from = estimate
  ) %>%
  extract(
    group_year,
    into = c("ethnicity", "socstatus", "year"),
    regex = "^(?:asian|black|white|mixed)_(.*)_(sc[1-7])_(\\d{4})$",
    remove = FALSE
  ) %>%
  mutate(
    group = paste(ethnicity, socstatus, sep = "_"),
    year = as.integer(year)
  ) %>%
  relocate(group, ethnicity, socstatus, year, .after = group_year)

coef_df %>% pivot_longer(
  cols = c(beta_E,beta_S),
  names_to = "dimension",
  values_to = "coefficient"
) %>% filter(socstatus %in%  c("sc1","sc7")) %>%
  ggplot(aes(x = year, y = coefficient,
             color = socstatus,
             linetype = dimension)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Linear Coefficient",
    title = "Differences in linear coefficients SC1-SC7",
    color = "Social class"
  ) +
  scale_linetype_manual(
    name = "Dimension",
    values = c(
      "beta_E" = "solid",
      "beta_S" = "dashed"
    ),
    labels = c(
      "beta_E" = "Ethnic",
      "beta_S" = "Social class"
    )
  ) +
facet_wrap(~ ethnicity, scales = "free_y") +
  scale_x_continuous(breaks = c(2001, 2011, 2021)) +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  ) +
  guides(
    color = guide_legend(ncol = 2, order = 1),
    linetype = guide_legend(ncol = 2, order = 1)
  )
ggsave("sc1_sc7_linearcomp.jpg", width = 8, height = 6)

# Wald test for linear coefficient differences

# Function: Wald test for the difference in linear coefficients

test_linear_difference <- function(
    model,
    year,
    dimension = "frac_homoe",
    reference_grouplin,
    reference_groupoth,
    vcov_type = ~ borough
) {
  
  # 1. Extract coefficients and variance-covariance matrix
  cf <- coef(model)
  V <- vcov(model, vcov = vcov_type)
  
  # 2. Identify linear coefficient for Indian SC1
  name_ref <- paste0(
     reference_grouplin , # "group_year::asian_indian_sc1_",
    year, ":", dimension
  )
  
  # 3. Identify linear coefficient for Indian SC7
  name_other <- paste0(
    reference_groupoth , # "group_year::asian_indian_sc7_",
    year, ":", dimension
  )
  
  # 4. Select the two linear coefficients
  terms <- c(name_ref, name_other)
  
  # Check coefficient availability
  if (!all(terms %in% names(cf))) {
    stop(
      "Missing coefficients:\n",
      paste(
        setdiff(terms, names(cf)),
        collapse = "\n"
      )
    )
  }
  
  # 5. Extract coefficient estimates
  b <- cf[terms]
  
  # 6. Extract their covariance matrix
  V_sub <- V[terms, terms, drop = FALSE]
  
  # 7. Wald restriction:
  # beta_linear_SC1 - beta_linear_SC7 = 0
  R <- matrix(
    c(1, -1),
    nrow = 1
  )
  
  # 8. Calculate difference between linear coefficients
  difference <- as.numeric(
    R %*% b
  )
  
  # 9. Calculate variance and standard error
  variance <- as.numeric(
    R %*% V_sub %*% t(R)
  )
  
  se <- sqrt(variance)
  
  # 10. Calculate Wald statistic
  W <- (difference / se)^2
  
  # 11. Calculate p-value
  p <- pchisq(
    W,
    df = 1,
    lower.tail = FALSE
  )
  
  # 12. Return results
  tibble(
    year = year,
    dimension = dimension,
    
    beta_linear_ref = unname(b[1]),
    beta_linear_oth = unname(b[2]),
    
    diff_linear = difference,
    se_diff = se,
    
    ci_lower = difference - 1.96 * se,
    ci_upper = difference + 1.96 * se,
    
    Wald_chisq = W,
    df = 1,
    p_value = p,
    
    significant = p < 0.05
  )
}


# Run Wald tests for all years and dimensions, referencegroup and other

expand.grid(
  year = c(2001, 2011, 2021),
  dimension = c("frac_homoe", "frac_homos")
) %>%
  pmap_dfr(function(year, dimension) {
    
    test_linear_difference(
      model = mod_B_quad,
      year = year,
      reference_grouplin = "group_year::asian_indian_sc3_",
      reference_groupoth = "group_year::asian_indian_sc5_",
      dimension = dimension
    )
    
  })


#### Marginal effects quadratic model

coef_df_sq <- tibble(
  term = names(coef(mod_B_quad)),
  estimate = unname(coef(mod_B_quad))
) %>%
  extract(
    term,
    into = c("group_year", "dimension", "order"),
    regex = "^group_year::(.+):(frac_homoe|frac_homos)(_sq)?$",
    remove = FALSE
  ) %>%
  filter(!is.na(group_year)) %>%
  mutate(
    order = if_else(order == "_sq", "quadratic", "linear"),
    dimension = recode(
      dimension,
      frac_homoe = "Ethnic",
      frac_homos = "Social"
    )
  ) %>%
  select(group_year, dimension, order, estimate) %>%
  pivot_wider(
    names_from = order,
    values_from = estimate
  ) %>%
  extract(
    group_year,
    into = c("ethnicity", "socstatus", "year"),
    regex = "^(?:asian|black|white|mixed)_(.*)_(sc[1-7])_(\\d{4})$",
    remove = FALSE
  ) %>%
  mutate(
    group = paste(ethnicity, socstatus, sep = "_"),
    year = as.integer(year)
  ) %>%
  relocate(
    group, ethnicity, socstatus, year,
    .after = group_year
  ) 

# Variance-covariance matrix
V <- vcov(mod_B_quad, vcov = ~ borough)

marginal_effects <- coef_df_sq %>%
  crossing(concentration = seq(0, 100, 1)) %>%
  mutate(
    marginal_effect = linear + 2 * quadratic * concentration
  ) %>% 
  mutate(
    dimension_var = recode(
      dimension,
      Ethnic = "frac_homoe",
      Social = "frac_homos"
    ),
    
    inference = pmap(
      list(group_year, dimension_var, concentration),
      function(group_year, dimension_var, concentration) {
        
        # Identify linear and quadratic coefficients
        name_linear <- paste0(
          "group_year::", group_year,
          ":", dimension_var
        )
        
        name_quadratic <- paste0(
          "group_year::", group_year,
          ":", dimension_var, "_sq"
        )
        
        terms <- c(name_linear, name_quadratic)
        
        if (!all(terms %in% rownames(V))) {
          stop(
            "Missing coefficients: ",
            paste(setdiff(terms, rownames(V)),
                  collapse = ", ")
          )
        }
        
        # Extract 2x2 covariance matrix
        V_sub <- V[terms, terms, drop = FALSE]
        
        # Gradient of marginal effect
        R <- matrix(
          c(1, 2 * concentration),
          nrow = 1
        )
        
        # Variance and standard error
        variance <- as.numeric(
          R %*% V_sub %*% t(R)
        )
        
        tibble(
          se_margeff = sqrt(variance)
        )
      }
    )
  ) %>%
  unnest(inference) %>%
  mutate(
    z_value = marginal_effect / se_margeff,
    
    p_value_margeff = 2 * pnorm(-abs(z_value)),
    
    ci_lower_margeff = marginal_effect - 1.96 * se_margeff,
    ci_upper_margeff = marginal_effect + 1.96 * se_margeff,
    
    significant = p_value_margeff < 0.05
  ) %>%
  select(-dimension_var)

for (i in unique(coef_df_sq$ethnicity)) { 
  
marginal_effects %>%
  filter(
    ethnicity == i
  ) %>%
  ggplot(
    aes(
      x = concentration,
      y = marginal_effect,
      color = factor(year) ,
      linetype = dimension
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line() +
  facet_wrap(~ socstatus, scales = "free_y") +
  scale_linetype_manual(
    values = c(
      "Ethnic" = "solid",
      "Social" = "dashed"
    ),
    labels = c(
      "Ethnic" = "Ethnic",
      "Social" = "Socio-occupational"
    )
  ) +
    scale_color_manual(
      values = c(
        "2001" = "salmon",
        "2011" = "green",
        "2021" = "blue"
      )
    ) +
  labs(
    x = "Residential concentration (%)",
    y = "Marginal effect on log-odds",
    color = "Dimension",
    linetype = "Year",
    title = str_to_title(i)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  ) +
    guides(
      color = guide_legend(ncol = 3, order = 1),
      linetype = guide_legend(ncol = 2, order = 1)
    )
  ggsave(paste0("figures_margeff/",i,".jpg"), width = 8, height = 8)
}

# within-class ethnic Wald of marginal effects

# Variance-covariance matrix clustered by borough
V <- vcov(mod_B_quad, vcov = ~ borough)

# List to store results
ethgrw <- "indian"
wald_list <- list()

for (i in paste0("sc", 1:7)) {
  
  # Reference class
  reference_sc <- marginal_effects %>%
    filter(
      ethnicity == ethgrw , # "indian",
      socstatus == i
    ) %>%
    select(
      year, dimension, concentration,
      marginal_reference = marginal_effect
    )
  
  # Compare reference class with every other Indian class
  marginal_diff <- marginal_effects %>%
    filter(
      ethnicity == ethgrw, #"indian",
      socstatus != i
    ) %>%
    left_join(
      reference_sc,
      by = c("year", "dimension", "concentration")
    ) %>%
    mutate(
      reference_class = i,
      difference = marginal_reference  - marginal_effect
    )
  
  # Wald tests
  wald_marginal <- marginal_diff %>%
    mutate(
      dimension_var = recode(
        dimension,
        Ethnic = "frac_homoe",
        Social = "frac_homos"
      )
    ) %>%
    mutate(
      test = pmap(
        list(
          year,
          dimension_var,
          concentration,
          reference_class,
          socstatus,
          difference
        ),
        
        function(year, dimension_var, concentration,
                 reference_class, socstatus, difference) {
          
          # Reference coefficient names
          name_reference <- paste0(
           "group_year::",macreth,"_", ethgrw, "_",          #  "group_year::asian_indian_",
            reference_class, "_",
            year, ":", dimension_var
          )
          
          # Comparison coefficient names
          name_other <- paste0(
            "group_year::",macreth,"_", ethgrw, "_",          #   "group_year::asian_indian_",
            socstatus, "_",
            year, ":", dimension_var
          )
          
          # Four coefficients
          terms <- c(
            name_reference,
            paste0(name_reference, "_sq"),
            name_other,
            paste0(name_other, "_sq")
          )
          
          # Check coefficients
          if (!all(terms %in% rownames(V))) {
            stop(
              "Missing coefficients: ",
              paste(
                setdiff(terms, rownames(V)),
                collapse = ", "
              )
            )
          }
          
          # Covariance matrix
          V_sub <- V[terms, terms, drop = FALSE]
          
          # Wald restriction:
          # ME_reference(C) - ME_other(C)
          R <- matrix(
            c(
              1, 2 * concentration,
              -1, -2 * concentration
            ),
            nrow = 1
          )
          
          # Variance and standard error
          variance <- as.numeric(
            R %*% V_sub %*% t(R)
          )
          
          se <- sqrt(variance)
          
          # Wald statistic
          W <- (difference / se)^2
          
          tibble(
            se = se,
            ci_lower = difference - 1.96 * se,
            ci_upper = difference + 1.96 * se,
            Wald_chisq = W,
            p_value = pchisq(
              W,
              df = 1,
              lower.tail = FALSE
            )
          )
        }
      )
    ) %>%
    unnest(test) %>%
    select(-dimension_var)
  
  # Store results for this reference class
  wald_list[[i]] <- wald_marginal
}

# Combine all results
wald_marginal_all <- bind_rows(wald_list)
write_xlsx(wald_marginal_all, "wald_marginal_indian.xlsx")


#within-class ethnic Wald for marginal effects differences

# Variance-covariance matrix
V <- vcov(mod_B_quad, vcov = ~ borough)

# Keep only necessary columns
marginal_base <- marginal_effects %>%
  select(
    group_year,
    ethnicity,
    socstatus,
    year,
    dimension,
    concentration,
    marginal_effect
  )

# Create all unique ethnic-subgroup pairs
# within the same social class
marginal_diff_ethnic <- marginal_base %>%
  inner_join(
    marginal_base,
    by = c(
      "socstatus",  # the socstatus allow for the same social status 
      "year",
      "dimension",
      "concentration"
    ),
    suffix = c("_ref", "_other"),
    relationship = "many-to-many"
  ) %>%
  filter(
    group_year_ref < group_year_other # alphabetical order
  ) %>%
  mutate(
    difference = marginal_effect_ref -
      marginal_effect_other
  )

# Wald tests
wald_ethnic <- marginal_diff_ethnic %>%
  mutate(
    dimension_var = recode(
      dimension,
      Ethnic = "frac_homoe",
      Social = "frac_homos"
    )
  ) %>%
  mutate(
    test = pmap(
      list(
        group_year_ref,
        group_year_other,
        dimension_var,
        concentration,
        difference
      ),
      
      function(group_year_ref, group_year_other,
               dimension_var, concentration, difference) {
        
        # Reference coefficient
        name_reference <- paste0(
          "group_year::",
          group_year_ref,
          ":",
          dimension_var
        )
        
        # Comparison coefficient
        name_other <- paste0(
          "group_year::",
          group_year_other,
          ":",
          dimension_var
        )
        
        # Four coefficients
        terms <- c(
          name_reference,
          paste0(name_reference, "_sq"),
          name_other,
          paste0(name_other, "_sq")
        )
        
        # Verify coefficient availability
        if (!all(terms %in% rownames(V))) {
          stop(
            "Missing coefficients: ",
            paste(
              setdiff(terms, rownames(V)),
              collapse = ", "
            )
          )
        }
        
        # Covariance matrix
        V_sub <- V[terms, terms, drop = FALSE]
        
        # Wald restriction
        R <- matrix(
          c(
            1, 2 * concentration,
            -1, -2 * concentration
          ),
          nrow = 1
        )
        
        # Variance of marginal-effect difference
        variance <- as.numeric(
          R %*% V_sub %*% t(R)
        )
        
        se <- sqrt(variance)
        
        # Wald statistic
        W <- (difference / se)^2
        
        tibble(
          se_diff = se,
          ci_lower_diff = difference - 1.96 * se,
          ci_upper_diff = difference + 1.96 * se,
          Wald_chisq = W,
          p_value_diff = pchisq(
            W,
            df = 1,
            lower.tail = FALSE
          )
        )
      }
    )
  ) %>%
  unnest(test) %>%
  select(-dimension_var)

# ODDS RATIO of quadratic model

or_quadratic <- coef_df_sq %>%
  crossing(
    concentration = seq(0, 99, 1)
  ) %>%
  mutate(
    OR = exp(
      linear + quadratic * (2 * concentration + 1)
    ),
    pct_change_odds = (OR - 1) * 100
  )

# simulated probabilities for quadratic model

simulate_probabilities_quad <- function(
    model,
    data,
    selected_group,
    x_grid = seq(0, 100, 1)
) {

  # Extract linear and quadratic coefficients
  beta_E1 <- coef(model)[
    paste0("group_year::", selected_group, ":frac_homoe")
  ]

  beta_E2 <- coef(model)[
    paste0("group_year::", selected_group, ":frac_homoe_sq")
  ]

  beta_S1 <- coef(model)[
    paste0("group_year::", selected_group, ":frac_homos")
  ]

  beta_S2 <- coef(model)[
    paste0("group_year::", selected_group, ":frac_homos_sq")
  ]

  stopifnot(
    length(beta_E1) == 1, !is.na(beta_E1),
    length(beta_E2) == 1, !is.na(beta_E2),
    length(beta_S1) == 1, !is.na(beta_S1),
    length(beta_S2) == 1, !is.na(beta_S2)
  )

  # Observed choice set
  df_g <- data %>%
    filter(group_year == selected_group)

  df_g$pred_count <- as.numeric(
    predict(model, newdata = df_g, type = "response")
  )

  # Simulate concentrations
  map_dfr(x_grid, function(x) {

    probabilities <- sapply(seq_len(nrow(df_g)), function(j) {

      # ETHNIC CONCENTRATION
      mu_E <- df_g$pred_count

      delta_E <-
        beta_E1 * (x - df_g$frac_homoe[j]) +
        beta_E2 * (x^2 - df_g$frac_homoe_sq[j])

      # here the counterfactual multinomial probabilities for ethnic dimension
      mu_E[j] <- mu_E[j] * exp(delta_E)

      prob_E <- mu_E[j] / sum(mu_E)

      # SOCIO-OCCUPATIONAL CONCENTRATION
      mu_S <- df_g$pred_count

      delta_S <-
        beta_S1 * (x - df_g$frac_homos[j]) +
        beta_S2 * (x^2 - df_g$frac_homos_sq[j])

      # here the counterfactual multinomial probabilities for social class dimension
      mu_S[j] <- mu_S[j] * exp(delta_S)

      prob_S <- mu_S[j] / sum(mu_S)

      c(
        ethnic = prob_E,
        social = prob_S
      )
    })

    tibble(
      group_year = selected_group,
      concentration = x,
      # the average probability from boroughs is extracted
      mean_prob_ethnic = mean(probabilities["ethnic", ]),
      mean_prob_social = mean(probabilities["social", ])
    )
  })
}

dir.create("figures_prob_quad", showWarnings = FALSE)

groups <- unique(
  str_remove(df_choice$group_year, "_\\d{4}$")
)

for (i in groups) {

  selected_groups <- paste0(
    i, "_", c(2001, 2011, 2021)
  )

  sim_results_all <- map_dfr(
    selected_groups,
    ~ simulate_probabilities_quad(
      model = mod_B_quad,
      data = df_choice,
      selected_group = .x
    )
  )

  sim_results_long <- sim_results_all %>%
    mutate(
      year = as.integer(
        str_extract(group_year, "\\d{4}$")
      )
    ) %>%
    pivot_longer(
      cols = c(mean_prob_ethnic, mean_prob_social),
      names_to = "dimension",
      values_to = "probability"
    ) %>%
    mutate(
      dimension = recode(
        dimension,
        mean_prob_ethnic = "Ethnic concentration",
        mean_prob_social = "Social-class concentration"
      )
    )

  group_label <- i %>%
    str_remove("^(asian|black|white|mixed)_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Sc([1-7])$", "SC\\1")

  p <- ggplot(
    sim_results_long,
    aes(
      x = concentration,
      y = probability,
      color = factor(year),
      linetype = dimension
    )
  ) +
    geom_line(linewidth = 1.1) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_linetype_manual(
      name = "Dimension",
      values = c(
        "Ethnic concentration" = "solid",
        "Social-class concentration" = "dashed"
      ),
      labels = c(
        "Ethnic concentration" = "Ethnic",
        "Social-class concentration" = "Socio-occupational"
      )
    ) +
    labs(
      x = "Concentration (%)",
      y = "Average predicted residential probability",
      color = "Year",
      title = paste0(group_label, ": 2001–2021")
    ) +
    theme_bw()

  ggsave(
    filename = paste0(
      "figures_prob_quad/", i, ".jpg"
    ),
    plot = p,
    width = 9,
    height = 6,
    dpi = 300
  )
}