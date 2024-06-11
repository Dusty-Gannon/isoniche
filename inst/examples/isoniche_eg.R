
library(isoniche)
library(dplyr)
library(ggplot2)
set.seed(9986)

# Load the example data
data("eg_data")

# ggplot(dat, aes(x = V1, y = V2, color = grp)) +
#   geom_point()

# fitting the model, one bivariate normal
# distribution per group
mfit <- isoniche(
  mean = list(
    y1 ~ grp + x,
    y2 ~ grp + x
  ),
  var = list(
    ~ grp,
    ~ grp,
    ~ grp
  ),
  data = eg_data
)

# make new data for constructing ellipses
newdat <- data.frame(
  grp = unique(eg_data$grp),
  x = rep(0, length(unique(eg_data$grp)))
)

plotdf <- construct_ellipses(mfit, newdat, n = 50)

ggplot(data = eg_data, aes(x = y1, y = y2, color = grp)) +
  geom_point() +
  geom_path(
    data = plotdf,
    aes(group = ellipse_id),
    size = 0.2
  ) +
  theme_classic() +
  scale_color_manual(values = c("steelblue", "brown"))

sea_post <- sea(mfit, newdat, n = 250)

ggplot(sea_post, aes(x = sea)) +
  geom_density(aes(fill = grp), alpha = 0.5, color = "grey") +
  theme_classic() +
  scale_fill_manual(values = c("steelblue", "brown"))

