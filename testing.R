
library(rstan)
library(tidyverse)

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

ggplot(data = dat, aes(x = V1, y = V2, color = grp)) +
  geom_point(alpha = 0.5)


# construct data list
X_1 <- model.matrix(~ grp, data = dat)

grps <- length(unique(dat$grp))

datlist <- list(
  N = nrow(dat),
  P = c(grps, grps),
  K = c(grps, grps),
  J = grps,
  X = cbind(X_1, X_1),
  Z = cbind(X_1, X_1),
  G = X_1,
  y = as.matrix(select(dat, V1:V2))
)

stanmod <- isoniche(
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
test <- sampling(stanmod, data = datlist)

# construct standard ellipses
X_new <- unique(datlist$X)
Z_new <- unique(datlist$Z)
G_new <- unique(datlist$G)

summary(test, pars = "gamma")

B <- matrix(data = 0, nrow = sum(datlist$P), ncol = 2)
B[1:datlist$P[1], 1] <- summary(test, pars = "beta_1")$summary[, "mean"]
B[(datlist$P[1] + 1):nrow(B), 2] <- summary(test, pars = "beta_2")$summary[, "mean"]

Z <- matrix(data = 0, nrow = sum(datlist$K), ncol = 2)
Z[1:datlist$K[1], 1] <- summary(test, pars = "zeta_1")$summary[, "mean"]
Z[(datlist$K[1] + 1):nrow(Z), 2] <- summary(test, pars = "zeta_2")$summary[, "mean"]

gamma <- summary(test, pars = "gamma")$summary[, "mean"]

mu_new <- map(
  1:nrow(X_new),
  ~ {
    as.double(X_new[.x, ] %*% B)
  }
)

Omega_new <- map(
  1:nrow(G_new),
  ~ {
    Omega <- diag(nrow = 2)
    rho <- 2 * plogis(as.double(G_new[.x, ] %*% gamma)) - 1
    Omega[1, 2] <- rho
    Omega[2, 1] <- rho
    Omega
  }
)

sigma_new <- map(
  1:nrow(Z_new),
  ~ {
    exp(as.double(Z_new[.x, ] %*% Z))
  }
)

Sigma_new <- map2(
  sigma_new, Omega_new,
  ~ {diag(.x) %*% .y %*% diag(.x)}
)

L_new <- map(
  Sigma_new,
  ~ t(chol(.x))
)

theta <- rbind(
  cos(seq(0, 2 * pi, length.out = 100)),
  sin(seq(0, 2 * pi, length.out = 100))
)

grp1 <- matrix(nrow = 100, ncol = 2)
grp2 <- grp1

for(i in 1:100){
  grp1[i, ] <- t(mu_new[[1]] + L_new[[1]] %*% theta[, i])
  grp2[i, ] <- t(mu_new[[2]] + L_new[[2]] %*% theta[, i])
}

dat_ellipse <- as_tibble(rbind(grp1, grp2))
dat_ellipse <- dat_ellipse %>% mutate(
  grp = factor(rep(c(1, 2), each = 100))
)

scatter +
  geom_path(
    data = dat_ellipse,
    aes(x = V1, y = V2, color = grp)
  )
ggsave("testplot.png", height = 4, width = 5, units = "in")


