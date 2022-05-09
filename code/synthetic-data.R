library(lme4)
library(arm)
library(here)

## Multi-level power calculations from Gelman's ARM book:
BHET.fake <- function (mu_b3, mu_g3, mu_g4){
  wave <- rep (seq(1,3,1), 19*50) 
  person <- rep (1:19, each=3)           # person ID’s
  village <- rep (1:50, each=19*3)
  treatment <- sample (rep (0:1, 50/2))
  treated <- treatment[village]
  post <- ifelse(wave > 1, 1, 0)
  policy <- ifelse(post==1 & treated==1, 1, 0)

  # # hyperparameters:
  mu.b0.true <- 0        # more generally, these could
  sigma.b0.true <- 0.1
  
  mu.b1.true <- 0.25       # be specified as additional
  mu.b2.true <- 0.25       # arguments to the function
  mu.b3.true <- mu_b3
  
  mu.g0.true <- 0
  sigma.g0.true <- 0.1
  mu.g1.true <- 0.25
  mu.g2.true <- 0.25
  mu.g3.true <- mu_g3
  mu.g4.true <- mu_g4
  
  sigma.m.true <- 1
  sigma.y.true <- 1

  # # village-level parameters
  b0.true <- rnorm(50, mu.b0.true, sigma.b0.true) 
  b1.true <- rnorm(50, mu.b1.true, 0)
  b2.true <- rnorm(50, mu.b2.true, 0)
  b3.true <- rnorm(50, mu.b3.true, 0)
  
  g0.true <- rnorm(50, mu.g0.true, sigma.g0.true) 
  g1.true <- rnorm(50, mu.g1.true, 0)
  g2.true <- rnorm(50, mu.g2.true, 0)
  g3.true <- rnorm(50, mu.g3.true, 0)
  g4.true <- rnorm(50, mu.g4.true, 0)

  # # data
  mm <- b0.true[village] + b1.true[village]*treated
     + b2.true[village]*post + b3.true[village]*policy
  s <- sigma.m.true
  location <- log(mm^2 / sqrt(s^2 + mm^2))
  shape <- sqrt(log(1 + (s^2 / mm^2)))

  m <- rnorm (50*19*3, b0.true[village] + b1.true[village]*treated
    + b2.true[village]*post + b3.true[village]*policy, sigma.m.true)
  
  # m <- rlnorm(50*19*3, location, shape)
  
  y <- rnorm (50*19*3, g0.true[village] + g1.true[village]*treated 
    + g2.true[village]*post + g3.true[village]*policy 
    + g4.true[village]*m, sigma.y.true)
  
  return (data.frame (y, m, wave, person, village, treated, post, policy))
}

# generate fake data
bhet.fake <- BHET.fake(mu_b3=1, mu_g3=0.5, mu_g4=0.5)

# write to data folder
write.csv(bhet.fake, 
          file = here( "data-clean", "bhet-fake.csv"))
