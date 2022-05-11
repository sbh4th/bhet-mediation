library(tidyverse)
library(fixest)
library(modelsummary)
library(osfr)
library(brms)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(here)

# read in data
d <- read_csv(here("data-clean", "bhet-fake.csv"))

# DD model
m1 <- feols(y ~ policy | village + wave, data=d, cluster="village")

## Mediated model
m2 <- feols(y ~ policy + m | village + wave, data=d, cluster="village")

# Summary
modelsummary(list(m1, m2),
  gof_omit = 'DF|Deviance|R2|AIC|BIC|Log.Lik|ICC|RMSE')

# Bayesian setup
# First the outcome model
d1 <- d %>%
  mutate(g = as.factor(village), # dummies for village
         t = as.factor(wave))    # dummies for wave

y_model <- bf(y ~ m + policy + g + t)

# now the mediator model
m_model <- bf(m ~ policy + g + t )

# multivariate model
bm <-
  brm(data = d1, 
      family = gaussian,
      y_model + m_model + set_rescor(FALSE),
      cores = 4, backend = "cmdstanr")

# pull posterior samples
post <- posterior_samples(bm)

# Define and plot effects ----

## define counterfactual quantities ====
post <-
  post %>% 
  mutate(cde = b_y_policy,
         nie = b_y_m * b_m_policy,
         te = cde + nie,
         pm = nie / te,
         pe = (te - cde) / te)



## plot of total effect and cde ====

post %>% 
  ggplot(aes(x = te)) + 
    geom_density(color="#1b9e77", fill="#1b9e77", alpha=0.2) + 
    geom_density(aes(x=cde), 
                 color="#d95f02", fill="#d95f02", alpha=0.2) +
    xlab("Hypothetical effect size") + ylab("Posterior density") + 
    annotate("text", x = 1, y = 4.5, 
           label="Total effect", size = 6, color = '#1b9e77') +
    annotate("text", x = 0.48, y = 5, 
           label="CDE", size = 6, color = '#d95f02') +
    theme_classic() + theme(axis.title = element_text(size=18))


## table of estimates ====
post %>%
  select(cde, nie, te, pm, pe) %>%
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>%
  median_qi()


