
library(rstan)
library(tidyverse)
library(isoniche)
set.seed(9986)

# simulating some example data
rho <- c(0, 0.5)
Omega_1 <- rbind(
  c(1, rho[1]),
  c(rho[1], 1)
)
Omega_2 <- rbind(
  c(1, rho[2]),
  c(rho[2], 1)
)
sigma_1 <- c(0.5, 0.5)
sigma_2 <- c(0.2, 1)

Sigma_1 <- diag(sigma_1) %*% Omega_1 %*% diag(sigma_1)
Sigma_2 <- diag(sigma_2) %*% Omega_2 %*% diag(sigma_2)

y_1 <- mvtnorm::rmvnorm(100, sigma = Sigma_1)
y_2 <- mvtnorm::rmvnorm(100, mean = c(1, 2), sigma = Sigma_2)

# construct the df
dat <- as_tibble(rbind(y_1, y_2))
dat <- dat %>% mutate(
  grp = factor(rep(c(1,2), each = 100))
)

# fitting the model, one bivariate normal
# distribution per group
mfit <- isoniche(
  mean = list(
    V1 ~ grp,
    V2 ~ grp
  ),
  var = list(
    ~ grp,
    ~ grp,
    ~ grp
  ),
  data = dat,
  cores = 4
)

# make new data for constructing ellipses
newdat <- data.frame(
  grp = unique(dat$grp)
)

plotdf <- construct_ellipses(mfit, newdat)

ggplot(data = dat, aes(x = V1, y = V2, color = grp)) +
  geom_point() +
  geom_path(
    data = plotdf,
    aes(group = ellipse_id)
  )


