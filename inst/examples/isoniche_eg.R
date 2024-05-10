
library(isoniche)
# library(dplyr)
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

# Continuous variable influencing the mean
x <- rnorm(200, sd = 2)
mu_1 <- diag(c(0.1, 0.8)) %*% rbind(x[1:100], x[1:100])
mu_2 <- diag(c(0.1, 0.8)) %*% rbind(x[101:200], x[101:200])

y_1 <- sapply(
  1:100,
  function(i, mu, sigma){
    mvtnorm::rmvnorm(1, as.double(mu[, i]), sigma = sigma)
  },
  mu = mu_1, sigma = Sigma_1
) %>% t()

y_2 <- sapply(
  1:100,
  function(i, mu, sigma){
    mvtnorm::rmvnorm(1, c(1, 2) + as.double(mu[, i]), sigma = sigma)
  },
  mu = mu_2, sigma = Sigma_2
) %>% t()

# construct the df
dat <- as_tibble(rbind(y_1, y_2))
dat <- dat %>% mutate(
  grp = factor(rep(c(1,2), each = 100)),
  x = x
)

# ggplot(dat, aes(x = V1, y = V2, color = grp)) +
#   geom_point()

# fitting the model, one bivariate normal
# distribution per group
mfit <- isoniche(
  mean = list(
    V1 ~ grp + x,
    V2 ~ grp + x
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
  grp = unique(dat$grp),
  x = rep(0, length(unique(dat$grp)))
)

plotdf <- construct_ellipses(mfit, newdat, n = 10)

ggplot(data = dat, aes(x = V1, y = V2, color = grp)) +
  geom_point() +
  geom_path(
    data = plotdf,
    aes(group = ellipse_id),
    size = 0.2
  ) +
  theme_classic() +
  scale_color_manual(values = c("steelblue", "brown"))


