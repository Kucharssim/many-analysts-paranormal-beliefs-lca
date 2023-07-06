library(here)
library(rstan)
library(tidyverse)
library(ggthemes)
source(here::here("R", "helpers.R"))
options(mc.cores = 8L) # for mclapply
rstan::rstan_options(auto_write = TRUE)
ggplot2::theme_set(theme_classic(base_size = 15))
plots <- list()
# data wrangling ----
source("item-dictionary.R")
df <- read.csv(here::here("Data", "data_paranormal_mma.csv"))
df$ppt <- seq_len(nrow(df))

if(interactive()) {
  View(df[is.na(df$weight),])
  View(df[!complete.cases(dplyr::select(df, item_dictionary$item)),])
}

df <- dplyr::select(df, ppt, weight, item_dictionary$item)

df$weight[is.na(df$weight)] <- 1

nrow(df)
sum(df$weight)


df_long <- tidyr::pivot_longer(df, cols = item_dictionary$item, names_to = "item", values_to = "response") |>
  dplyr::select(ppt, item, response, weight) |>
  dplyr::filter(!is.na(response)) |>
  dplyr::mutate(
    item = factor(item, levels = item_dictionary$item),
    label = factor(item, levels = item_dictionary$item, labels = item_dictionary$label),
    response = factor(response, levels = 1:5, labels = c("Definitely not", "Probably not", "Maybe", "Probably", "Definitely"), ordered = TRUE)
    )


# descriptives ----
global_labeller <- ggplot2::labeller(
  item = setNames(item_dictionary$label, item_dictionary$item)
)


# barplots of observed responses
df_long |>
  dplyr::group_by(item, response) |>
  dplyr::summarise(count = length(weight), weighted = sum(weight)) |>
  dplyr::ungroup() |>
  tidyr::pivot_longer(cols = c("count", "weighted"), names_to = "Type") |>
  ggplot2::ggplot(ggplot2::aes(x = response, y = value, group = Type, fill = Type)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  ggplot2::facet_wrap(~item, ncol = 6, labeller = global_labeller) +
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
  ggplot2::xlab("Response") +
  ggplot2::ylab("Counts") +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 10)) +
  ggthemes::scale_fill_colorblind()
# not much going on there, except that there are many responses "definitely not"
# -> this is consistent with a group of "sceptics"

df_long |>
  dplyr::group_by(ppt) |>
  dplyr::summarise(`Definitely not` = sum(response == "Definitely not"),
                   `Probably not` = sum(response %in% c("Definitely not", "Probably not"))) |>
  dplyr::ungroup() |>
  tidyr::pivot_longer(cols = c("Definitely not", "Probably not"), names_to = "Response", names_ptypes = list(Response = factor(ordered = TRUE))) |>
  ggplot2::ggplot(ggplot2::aes(x = value, group = Response, fill = Response)) +
  ggplot2::geom_histogram(position = ggplot2::position_dodge(width = 0.5), bins = 24) +
  ggthemes::scale_fill_colorblind() +
  ggplot2::xlab("Number of responses") +
  ggplot2::ylab("Counts") +
  ggplot2::ggtitle("Histogram of responses smaller than the given category per respondent")

# There is a large group of respondents who answered most of the items "Definitely not" and below.
# There appears a little bump of respondents who answer about a half of questions below "Definitely not"
# which could indicate some kind of clustering of item types... However, it is far from being clear to claim
# anything like this - it could also be just noise or an artifact.


# fitting models ----
## uninformed by CFA results ----
### constrained ----
stan_model <- rstan::stan_model(here::here("stan", "constrained.stan"))

stan_data  <- list(
  N_obs = nrow(df_long),
  N_itm = nrow(item_dictionary),
  N_ppt = length(unique(df_long$ppt)),
  N_cat = 5,
  y     = df_long$response |> as.integer(),
  ppt   = df_long$ppt |> as.factor() |> as.integer(),
  itm   = df_long$item |> as.factor() |> as.integer(),
  weight = (df_long |> group_by(ppt) |> summarise(weight = mean(weight)) |> dplyr::select(weight))$weight
)


for(class in 1:2) {
  stan_data$N_cls <- class
  stan_fit <- repeatOptim(
    n = 25,
    parallel = TRUE,
    object = stan_model,
    data = stan_data,
    draws = 1000,
    importance_resampling = TRUE,
    as_vector = FALSE
  )
  saveRDS(stan_fit, here::here("models", sprintf("constrained-%s.Rds", class)))
}
stan_data$N_cls <- 1

rstan::optimizing(  object = stan_model,
                    data = stan_data,
                    draws = 1000,
                    importance_resampling = TRUE,
                    as_vector = FALSE)
### unconstrained ----
stan_model <- stan_model(here::here("stan", "unconstrained.stan"))

stan_data  <- list(
  N_obs = nrow(df_long),
  N_itm = nrow(item_dictionary),
  N_ppt = length(unique(df_long$ppt)),
  N_cat = 5,
  y     = df_long$response |> as.integer(),
  ppt   = df_long$ppt |> as.factor() |> as.integer(),
  itm   = df_long$item |> as.factor() |> as.integer(),
  weight = (df_long |> group_by(ppt) |> summarise(weight = mean(weight)) |> dplyr::select(weight))$weight
)

for(N_cls in 1:16) {
  cat("Fitting", N_cls, "classes\n")
  stan_data$N_cls <- N_cls
  stan_fit <- repeatOptim(
    n = 25,
    parallel = TRUE,
    object = stan_model,
    data = stan_data,
    draws = 1000,
    importance_resampling = TRUE,
    as_vector = FALSE)
  saveRDS(stan_fit, here::here("models", sprintf("unconstrained-%s.Rds", N_cls)))
}

bic <- numeric()
for(class in 1:16) {
  stan_fit <- readRDS(here::here("models", sprintf("unconstrained-%s.Rds", class)))
  plot(stan_fit$refits$logLik)
  title(class)
  bic <- c(bic, stan_fit$par$bic)
}
# we could not get > 6 classes to converge

stan_eta <- data.frame(class = NULL, resp = NULL, estimate = NULL, lower = NULL, upper = NULL)
for (class in 1:N_cls) {
  for (resp in seq_len(nrow(item_dictionary))) {
    estimate <- stan_fit$par$eta[resp, class]
    draws <- stan_fit$theta_tilde[, sprintf("eta[%s,%s]", resp, class)]
    lower <- quantile(draws, 0.005)
    upper <- quantile(draws, 0.995)
    stan_eta <- rbind(stan_eta, data.frame(class, resp, estimate, lower, upper))
  }
}
unconstrained <- readRDS(here::here("models", "unconstrained-6.Rds"))

bic <- list()
files <- list.files("models")
for (file in files) {
  stan_fit <- readRDS(here::here("models", file))
  bic[[file]] <- stan_fit$par$bic
}
bic <- unlist(bic)

bic2 <- bic[
  c("cfa-constrained.Rds", "cfa-unconstrained.Rds",
    "constrained1.Rds", "constrained2.Rds",
    "unconstrained-1.Rds", "unconstrained-2.Rds")
]
bic2 <- bic2 - min(bic2)
(exp(-bic2/2) / sum(exp(-bic2/2))) |> barplot()


## informed by CFA results ----

### constrained ----
stan_model <- rstan::stan_model(here::here("stan", "cfa-constrained.stan"))

N_cls <- 2^4
stan_data <- list(
  N_obs = nrow(df_long),
  N_itm = nrow(item_dictionary),
  N_ppt = length(unique(df_long$ppt)),
  N_cat = 5,
  N_cls = N_cls,
  y     = df_long$response |> as.integer(),
  ppt   = df_long$ppt |> as.factor() |> as.integer(),
  itm   = df_long$item |> as.factor() |> as.integer(),
  weight = (df_long |> group_by(ppt) |> summarise(weight = mean(weight)) |> dplyr::select(weight))$weight,

  N_fac = 4,
  fac = dplyr::left_join(df_long, item_dictionary, by = "item")$factor |> as.integer(),
  N_lvl = 2
)
cls <- tidyr::expand_grid(eha=1:2, sr=1:2, ub=1:2, es=1:2)
lvl <- cls[,stan_data$fac] |> as.matrix() |> t()
stan_data$lvl <- lvl

stan_fit <- repeatOptim(
  n = 30,
  parallel = TRUE,
  object = stan_model,
  data = stan_data,
  draws = 1000,
  importance_resampling = TRUE,
  as_vector = FALSE
)
saveRDS(stan_fit, here::here("models", "cfa-constrained.Rds"))
rm(stan_model, stan_fit)

### semi-constrained  -----
stan_model <- rstan::stan_model(here::here("stan", "cfa-semi-constrained.stan"))

N_cls <- 2^4
stan_data <- list(
  N_obs = nrow(df_long),
  N_itm = nrow(item_dictionary),
  N_ppt = length(unique(df_long$ppt)),
  N_cat = 5,
  N_cls = N_cls,
  y     = df_long$response |> as.integer(),
  ppt   = df_long$ppt |> as.factor() |> as.integer(),
  itm   = df_long$item |> as.factor() |> as.integer(),
  weight = (df_long |> group_by(ppt) |> summarise(weight = mean(weight)) |> dplyr::select(weight))$weight,

  N_fac = 4,
  fac = dplyr::left_join(df_long, item_dictionary, by = "item")$factor |> as.integer(),
  N_lvl = 2
)
cls <- tidyr::expand_grid(eha=1:2, sr=1:2, ub=1:2, es=1:2)
lvl <- cls[,stan_data$fac] |> as.matrix() |> t()
stan_data$lvl <- lvl

stan_fit <- repeatOptim(
  n = 10,
  parallel = TRUE,
  object = stan_model,
  data = stan_data,
  draws = 1000,
  importance_resampling = TRUE,
  as_vector = FALSE
)
saveRDS(stan_fit, here::here("models", "cfa-semi-constrained.Rds"))
rm(stan_model, stan_fit)

### unconstrained ----
stan_model <- rstan::stan_model(here::here("stan", "cfa-unconstrained.stan"))

N_cls <- 2^4
stan_data <- list(
  N_obs = nrow(df_long),
  N_itm = nrow(item_dictionary),
  N_ppt = length(unique(df_long$ppt)),
  N_cat = 5,
  N_cls = N_cls,
  y     = df_long$response |> as.integer(),
  ppt   = df_long$ppt |> as.factor() |> as.integer(),
  itm   = df_long$item |> as.factor() |> as.integer(),
  weight = (df_long |> group_by(ppt) |> summarise(weight = mean(weight)) |> dplyr::select(weight))$weight,

  N_fac = 4,
  fac = dplyr::left_join(df_long, item_dictionary, by = "item")$factor |> as.integer(),
  N_lvl = 2
)
cls <- tidyr::expand_grid(eha=1:2, sr=1:2, ub=1:2, es=1:2)
lvl <- cls[,stan_data$fac] |> as.matrix() |> t()
stan_data$lvl <- lvl

stan_fit <- repeatOptim(
  n = 20,
  parallel = TRUE,
  object = stan_model,
  data = stan_data,
  draws = 1000,
  importance_resampling = TRUE,
  as_vector = FALSE
)
saveRDS(stan_fit, here::here("models", "cfa-unconstrained.Rds"))
rm(stan_model, stan_fit)

rstan::optimizing(  object = stan_model,
                    data = stan_data,
                    draws = 1000,
                    importance_resampling = TRUE,
                    as_vector = FALSE, verbose = TRUE)



# inference ----
## model comparison ----

models <- tibble::tibble(
  type = c(rep("Uninformed", 4), rep("Informed by CFA", 3)),
  name = c("Constrained(1)", "Constrained(2)", "Unconstrained(1)", "Unconstrained(2)", "Constrained", "Semi-constrained", "Unconstrained"),
  id = c("constrained-1", "constrained-2", "unconstrained-1", "unconstrained-2", "cfa-constrained", "cfa-semi-constrained", "cfa-unconstrained"),
  pi = c(1, 2, 1, 2, 16, 16, 16),
  eta = c(1, 2, 24, "$2 \\times 24$", 2, "$2 \\times 4$", "$2 \\times 24$"),
  c = "$4 \\times 24$",
  npar = c(73, 75, 96, 121, 89, 95, 135),
  logLik = NA,
  aic = NA,
  bic = NA
)

models_list <- lapply(models$id, function(id) {
  path <- here::here("models", paste0(id, ".Rds"))
  readRDS(path)
})

n_eff <- (df_long |> group_by(ppt) |> summarise(w = weight[1]))$w |> sum()
models$logLik <- sapply(models_list, "[[", "value")
models$aic <- with(models, 2 * npar - 2*logLik)
models$bic <- with(models, npar * log(n_eff) - 2 * logLik)

m_aic <- min(models$aic)
models$w_aic <- with(models, exp(-(aic-m_aic)/2) / sum(exp(-(aic-m_aic)/2)))
m_bic <- min(models$bic)
models$w_bic <- with(models, exp(-(bic-m_bic)/2) / sum(exp(-(bic-m_bic)/2)))

saveRDS(object = models, file = here::here("models", "model-comparison.Rds"))

plots[["model_comparison"]] <- models |>
  tidyr::pivot_longer(cols = c("aic", "bic"), values_to = "value", names_to = "crit") |>
  dplyr::mutate(crit = as.factor(crit), type = factor(type, levels = c("Uninformed", "Informed by CFA"))) |>
  ggplot2::ggplot(mapping = ggplot2::aes(x = name, y = value, group = crit, col = crit)) +
  ggplot2::geom_point(ggplot2::aes(y = value), size = 5) +
  ggplot2::facet_wrap(~type, scales = "free_x") +
  ggthemes::scale_color_colorblind(name = "", labels = c("AIC", "BIC")) +
  ggplot2::ylab(NULL) +
  ggplot2::xlab("Model") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(
  plot = plots[["model_comparison"]],
  filename = "model-comparison.eps",
  path = here::here("figures"),
  device = cairo_ps,
  width = 7, height = 5)

## cfa - unconstrained ----
lower_p <- 0.025
upper_p <- 0.975
cfa_unconstrained <- readRDS(here::here("models", "cfa-unconstrained.Rds"))
par <- cfa_unconstrained$par
par_samples <- cfa_unconstrained$theta_tilde

pars <- data.frame(
  par = colnames(cfa_unconstrained$theta_tilde),
  estimate = sapply(colnames(cfa_unconstrained$theta_tilde), function(p) {
    eval(parse(text = sprintf("par$%s", p)))
  }),
  lower = apply(cfa_unconstrained$theta_tilde, 2, quantile, lower_p),
  upper = apply(cfa_unconstrained$theta_tilde, 2, quantile, upper_p)
)

### eta ----
plots[["eta"]] <- pars |>
  dplyr::filter(startsWith(par, "eta")) |>
  tidyr::separate(par, into = c("par", "item", "class"), extra = "drop") |>
  dplyr::mutate(
    item = as.factor(as.integer(item)),
    class = factor(class, labels = c("Low", "High"))
    ) |>
  ggplot2::ggplot(ggplot2::aes(x = item, group = class)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, group = class, fill = class), alpha = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = estimate), linewidth = 1) +
  ggthemes::scale_fill_colorblind(name = "Level") +
  ggplot2::xlab("Item") +
  ggplot2::ylab(expression(eta))
ggplot2::ggsave(
  plot = plots[["eta"]],
  filename = "eta.eps",
  path = here::here("figures", "cfa-unconstrained"),
  device = cairo_ps,
  width = 7, height = 7)

### eta per factor ---
plots[["eta-grouped"]] <- pars |>
  dplyr::filter(startsWith(par, "eta")) |>
  tidyr::separate(par, into = c("par", "item", "class"), extra = "drop") |>
  dplyr::mutate(
    item = factor(as.integer(item), labels = sprintf("par_%s", unique(as.integer(item)))),
    class = factor(class, labels = c("Low", "High"))
  ) |>
  dplyr::left_join(item_dictionary, by = "item") |>
  dplyr::mutate(
    item = as.integer(item) |> factor(labels = 1:24)
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = item, group = class)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, group = class, fill = class), alpha = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = estimate), linewidth = 1) +
  ggthemes::scale_fill_colorblind(name = "Level") +
  ggplot2::xlab("Item") +
  ggplot2::ylab(expression(eta)) +
  ggplot2::facet_wrap(~factor, ncol = 4, scales = "free_x")
ggplot2::ggsave(
  plot = plots[["eta-grouped"]],
  filename = "eta-grouped.eps",
  path = here::here("figures", "cfa-unconstrained"),
  device = cairo_ps,
  width = 12, height = 8)

### cutpoints ----
plots[["cutpoints"]] <- pars |>
  dplyr::filter(startsWith(par, "cutpoints[")) |>
  tidyr::separate(par, into = c("par", "item", "index"), extra = "drop") |>
  dplyr::mutate(item = as.factor(as.integer(item))) |>
  ggplot2::ggplot(ggplot2::aes(x = item, group = index)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, group = index), alpha = 0.7) +
  ggplot2::geom_line(ggplot2::aes(y = estimate), linewidth = 1) +
  ggplot2::xlab("Item") +
  ggplot2::ylab("Cutpoints")
ggplot2::ggsave(
  plot = plots[["cutpoints"]],
  filename = "cutpoints.eps",
  path = here::here("figures", "cfa-unconstrained"),
  device = cairo_ps,
  width = 12, height = 7)

patchwork::wrap_plots(plots[c("eta", "cutpoints")], ncol = 1)
ggplot2::ggsave(
  filename = "eta-cutpoints.eps",
  path = here::here("figures", "cfa-unconstrained"),
  device = cairo_ps,
  width = 8, height = 8)

### probs ----
probs <- data.frame(item = NULL, response = NULL, level = NULL, estimate = NULL, lower = NULL, upper = NULL)
for(i in 1:24) {
  for(r in 1:5) {
    for(l in 1:2) {
      uc <- if(r == 1) {-Inf} else {par$cutpoints[i, r-1]}
      lc <- if(r == 5) { Inf} else {par$cutpoints[i, r]}
      et <- par$eta[i,l]
      estimate <- inv_logit(et - uc) - inv_logit(et - lc)

      uc <- if(r == 1) {-Inf} else {par_samples[, sprintf("cutpoints[%s,%s]", i, r-1)]}
      lc <- if(r == 5) { Inf} else {par_samples$theta_tilde[, sprintf("cutpoints[%s,%s]", i, r  )]}
      et <- par_samples[, sprintf("eta[%s,%s]", i, l)]
      est_ <- inv_logit(et - uc) - inv_logit(et - lc)
      lower <- quantile(est_, lower_p)
      upper <- quantile(est_, upper_p)
      probs <- rbind(
        probs,
        data.frame(
          item = i, response = r, level = l,
          estimate = estimate, lower = lower, upper = upper
        )
      )
    }
  }
}

plots[["probs"]] <- probs |>
  dplyr::mutate(
    level = factor(level, labels = c("Low", "High")),
    item = factor(item, labels = levels(item_dictionary$item)),
    response = ordered(response, labels = levels(df_long$response))
    ) |>
  ggplot2::ggplot(ggplot2::aes(x = response, y = estimate, group = level, fill = level)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.5) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~item, ncol = 6, labeller = global_labeller) +
  ggthemes::scale_fill_colorblind(name = "Level") +
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
  ggplot2::xlab("Response") +
  ggplot2::ylab("P(Response)")
ggplot2::ggsave(
  plot = plots[["probs"]],
  filename = "probs.eps",
  path = here::here("figures", "cfa-unconstrained"),
  device = cairo_ps,
  width = 14, height = 14
)

### marginal probs ----
marginal_probs <- data.frame(item = NULL, response = NULL, estimate = NULL, lower = NULL, upper = NULL)
for(i in 1:24) {
  for(r in 1:5) {
    est <- 0
    est_samples <- 0
    for(c in 1:16) {
      cp <- par$class_prob[c]
      uc <- if(r == 1) {-Inf} else {par$cutpoints[i, r-1]}
      lc <- if(r == 5) { Inf} else {par$cutpoints[i, r]}
      f <- item_dictionary[i, "factor"] |> as.integer()
      l <- cls[c,f] |> as.numeric()
      et <- par$eta[i,l]
      est <- est + cp * (inv_logit(et - uc) - inv_logit(et - lc))


      cp <- par_samples[, sprintf("class_prob[%s]", c)]
      uc <- if(r == 1) {-Inf} else {par_samples[, sprintf("cutpoints[%s,%s]", i, r-1)]}
      lc <- if(r == 5) { Inf} else {par_samples[, sprintf("cutpoints[%s,%s]", i, r  )]}
      et <- par_samples[, sprintf("eta[%s,%s]", i, l)]
      est_samples <- est_samples + cp * (inv_logit(et - uc) - inv_logit(et - lc))
    }
    marginal_probs <- rbind(
      marginal_probs,
      data.frame(
        item = i,
        response = r,
        estimate = est,
        lower = quantile(est_samples, lower_p),
        upper = quantile(est_samples, upper_p))
    )
  }
}

marginal_probs <- marginal_probs |>
  dplyr::mutate(
    item = factor(item, labels = levels(item_dictionary$item)),
    response = factor(response, labels = levels(df_long$response), ordered = TRUE)
  )

marginal_probs |>
  ggplot2::ggplot(ggplot2::aes(x = response)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper)) +
  ggplot2::geom_point(ggplot2::aes(y = estimate)) +
  #ggplot2::geom_bar(ggplot2::aes(y = estimate), stat = "identity") +
  ggplot2::facet_wrap(~item, ncol = 6, labeller = global_labeller) +
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 10)) +
  ggplot2::xlab("Response") +
  ggplot2::ylab("P(Response)")


plots[["marginal_probs"]] <- df_long |>
  dplyr::group_by(item, response) |>
  dplyr::summarise(observed_prop = length(weight), weighted_prop = sum(weight)) |>
  dplyr::ungroup() |>
  dplyr::group_by(item) |>
  dplyr::mutate(observed_prop = observed_prop / sum(observed_prop),
                weighted_prop = weighted_prop / sum(weighted_prop)) |>
  dplyr::ungroup() |>
  dplyr::left_join(marginal_probs, by = c("item", "response")) |>
  tidyr::pivot_longer(cols = c("observed_prop", "weighted_prop", "estimate"), names_to = "Type") |>
  dplyr::mutate(
    Type = factor(Type, levels = c("observed_prop", "weighted_prop", "estimate"), labels = c("Observed", "Weighted", "Estimated"))
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = response, y = value, group = Type, fill = Type)) +
  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
  ggplot2::facet_wrap(~item, ncol = 6, labeller = global_labeller) +
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
  ggplot2::xlab("Response") +
  ggplot2::ylab("P(Response)") +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 10)) +
  ggthemes::scale_fill_colorblind()
ggplot2::ggsave(
  plot = plots[["marginal_probs"]],
  filename = "marginal-probs.eps",
  path = here::here("figures", "cfa-unconstrained"),
  device = cairo_ps,
  width = 14, height = 14
)

### class sizes ----
classification <- apply(par$class_membership, 1, which.max) |> factor(levels = 1:16)
header <- c(2, 2, 1)
names(header) <- c(" ", sprintf("%s%% CI", (upper_p - lower_p) * 100), " ")
pars |>
  dplyr::filter(startsWith(par, "class_prob")) |>
  dplyr::mutate(
    par   = factor(par, levels = sprintf("class_prob[%s]", 1:16), labels = 1:16),
    label = c(
      "Skeptics",
      "ES",
      "UB",
      "ES + UB",
      "SR",
      "SR + ES",
      "SR + UB",
      "SR + UB + ES",
      "EHA",
      "EHA + ES",
      "EHA + UB",
      "EHA + UB + ES",
      "EHA + SR",
      "EHA + SR + ES",
      "EHA + SR + UB",
      "Believers"
    ),
    count = table(classification),
    prop = count / length(classification)
  ) |>
  dplyr::arrange(by = desc(estimate)) |>
  dplyr::select(label, estimate, lower, upper, count) |>
  knitr::kable(format = "html", digits = 2,
               row.names = FALSE, col.names = c("Class", "Proportion", "Lower", "Upper", "Count"),
               booktabs = TRUE,
               linesep = c(""),
               caption = "Proportions of the 16 classes estimated by the model. ") |>
  kableExtra::add_header_above(header = header) |>
  kableExtra::add_footnote("EHA - Belief in extraordinary human abilities; SR - Belief in supernatural reality; UB - Belief in unearthly beings; ES - Belief in everyday superstition.", notation = "symbol") |>
  kableExtra::add_footnote("Count - The number of participants classified using the maximum aposteriori principle. The observed counts are not exactly proportional to the estimated proportions because the model was weighted by the survey weights.", notation = "symbol")

# |>
#   ggplot2::ggplot(ggplot2::aes(x = par, y = estimate, ymin = lower, ymax = upper)) +
#   ggplot2::geom_errorbar(width = 0.2) +
#   ggplot2::geom_point() +
#   ggplot2::xlab("Class") +
#   ggplot2::ylab("Size")

### factor sizes ----
factor_sizes <- data.frame(fac=NULL, estimate=NULL, lower=NULL, upper=NULL)
for(fac in colnames(cls)) {
  cls_ind <- which(cls[[fac]] == 2)
  estimate <- sum(par$class_prob[cls_ind])
  nm <- sprintf("class_prob[%s]", cls_ind)
  lower <- quantile(rowSums(cfa_unconstrained$theta_tilde[, nm]), 0.025)
  upper <- quantile(rowSums(cfa_unconstrained$theta_tilde[, nm]), 0.975)
  factor_sizes <- rbind(
    factor_sizes,
    data.frame(fac = fac, estimate = estimate, lower = lower, upper = upper)
  )
}
factor_sizes$fac <- factor(factor_sizes$fac,
                           levels = c("eha", "sr", "ub", "es"),
                           labels = c("Extraordinary human abilities", "Supernatural reality", "Unearthly beings", "Everyday superstition")
                           )

header <- c(2, 2)
names(header) <- c(" ", sprintf("%s%% CI", (upper_p - lower_p) * 100))
knitr::kable(factor_sizes, format = "latex", digits = 2,
             row.names = FALSE, col.names = c("Factor", "Proportion", "Lower", "Upper"),
             booktabs = TRUE,
             linesep = c(""),
             caption = "Proportion of people who are not entirely sceptical towards one of the four groups of items. ") |>
  kableExtra::add_header_above(header = header)





p_eta <- ggplot(stan_eta, aes(x = as.factor(resp), y = estimate, group = class)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(class)), alpha = 0.9) +
  geom_line(col = "black") +
  xlab("item") +
  ylab("estimate") +
  ggtitle("Locations") +
  scale_fill_discrete(name = "class")

ggplot(stan_eta, mapping = aes(x = class, y = estimate)) +
  geom_line() +
  facet_wrap(~resp)

stan_cutpoints <- data.frame(item = NULL, cutpoint = NULL, estimate = NULL, lower = NULL, upper = NULL)
for (item in seq_len(nrow(item_dictionary))) {
  for (cut in 1:4) {
    estimate <- stan_fit$par$cutpoints[item, cut]
    draws <- stan_fit$theta_tilde[, sprintf("cutpoints[%s,%s]", item, cut)]
    lower <- quantile(draws, 0.005)
    upper <- quantile(draws, 0.995)
    stan_cutpoints <- rbind(stan_cutpoints, data.frame(item = item, cutpoint = cut, estimate = estimate, lower = lower, upper = upper))
  }
}

p_cut <- ggplot(stan_cutpoints, aes(x = as.factor(item), y = estimate, group = cutpoint)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  geom_line(col = "black") +
  xlab("item") +
  ylab("estimate") +
  ggtitle("Cutpoints")


patchwork::wrap_plots(p_eta, p_cut, ncol = 1)


stan_probs <- data.frame(class = NULL, resp = NULL, category = NULL, estimate = NULL)
for (class in 1:N_cls) {
  for (resp in 1:28) {
    for (category in 1:5) {
      estimate <- stan_fit$par$probs[resp, class, category]
      stan_probs <- rbind(stan_probs, data.frame(class, resp, category, estimate))
    }
  }
}

ggplot(stan_probs, mapping = aes(x = category, y = estimate, group = class, col = as.factor(class))) +
  geom_line() +
  facet_wrap(~resp)


ggplot(df_long, aes(x = response)) +
  geom_histogram() +
  facet_wrap(~item)



## exploratory results ----
# read all unconstrained models
unconstr_files <- list.files(here::here("models"), pattern = "unconstrained-", full.names = TRUE)
fits <- lapply(unconstr_files, readRDS)

# check if we managed to replicate the best fit across different optimalization runs
converged <- sapply(fits, function(f) {
  ll <- f$refits$logLik[f$refits$return_code == 0]
  if (length(ll) <= 1) return(FALSE)
  max_ll <- max(ll)
  replicated <- round(max_ll, 1) == round(ll, 1)
  if (sum(replicated) > 1) return(TRUE)
  return(FALSE)
})

fits <- fits[converged]
classes <- sapply(fits, function(f) f$par$class_prob |> length())
npar <- (classes-1) + 24 * classes + 24 * 3
logLik <- sapply(fits, "[[", "value")
model_comparison <- data.frame(
  classes = classes,
  npar    = npar,
  logLik  = logLik,
  aic     = 2 * npar - 2*logLik,
  bic     = npar * log(n_eff) - 2 * logLik
)
model_comparison

pars <- data.frame(
  model = numeric(), par = character(), estimate = numeric(), lower = numeric(), upper = numeric()
)
for(fit in fits) {
  pars_samples <- fit$theta_tilde
  par <- fit$par
  tmp <- data.frame(
    model = fit$par$class_prob |> length(),
    par = colnames(pars_samples),
    estimate = sapply(colnames(pars_samples), function(p) {
      eval(parse(text = sprintf("par$%s", p)))
    }),
    lower = apply(pars_samples, 2, quantile, lower_p),
    upper = apply(pars_samples, 2, quantile, upper_p)
  )
  pars <- rbind(pars, tmp)
}
rownames(pars) <- NULL

plots$exploratory$eta_overview <- pars |>
  dplyr::filter(startsWith(par, "eta")) |>
  tidyr::separate(col = "par", into = c("par", "item", "class", NA)) |>
  dplyr::mutate(
    class = as.factor(class),
    item = as.integer(item),
    model = sprintf("Model %s", model)) |>
  ggplot2::ggplot(ggplot2::aes(x = item, y = estimate, group = class, col = class, fill = class)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~model) +
  ggplot2::xlab("Item") +
  ggplot2::ylab(expression(eta))
ggplot2::ggsave(
  plot = plots$exploratory$eta_overview,
  filename = "eta-overview.eps",
  path = here::here("figures", "exploratory"),
  device = cairo_ps,
  width = 10, height = 8)

### Model 6 ----
plots$exploratory[["eta-grouped-6"]] <- pars |>
  dplyr::filter(startsWith(par, "eta"), model == 6) |>
  tidyr::separate(par, into = c("par", "item", "class"), extra = "drop") |>
  dplyr::mutate(
    item = factor(as.integer(item), labels = sprintf("par_%s", unique(as.integer(item)))),
    class = as.factor(class)
  ) |>
  dplyr::left_join(item_dictionary, by = "item") |>
  dplyr::mutate(
    item = as.integer(item) |> factor(labels = 1:24)
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = item, group = class)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, group = class, fill = class), alpha = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = estimate), linewidth = 1) +
  ggthemes::scale_fill_colorblind(name = "Class") +
  ggplot2::xlab("Item") +
  ggplot2::ylab(expression(eta)) +
  ggplot2::facet_wrap(~factor, ncol = 4, scales = "free_x")
ggplot2::ggsave(
  plot = plots$exploratory[["eta-grouped-6"]],
  filename = "eta-grouped-6.eps",
  path = here::here("figures", "exploratory"),
  device = cairo_ps,
  width = 11, height = 6)

plots$exploratory[["eta-averaged-6"]] <- pars |>
  dplyr::filter(startsWith(par, "eta"), model == 6) |>
  tidyr::separate(par, into = c("par", "item", "class"), extra = "drop") |>
  dplyr::mutate(
    item = factor(as.integer(item), labels = sprintf("par_%s", unique(as.integer(item)))),
    class = as.factor(class)
  ) |>
  dplyr::left_join(item_dictionary, by = "item") |>
  dplyr::mutate(
    item = as.integer(item) |> factor(labels = 1:24)
  ) |>
  dplyr::group_by(factor, class) |>
  dplyr::summarise(mean = mean(estimate)) |>
  ggplot2::ggplot(ggplot2::aes(x = factor, y = mean, group = class, col = class)) +
  ggplot2::geom_line(linewidth = 1) +
  ggthemes::scale_color_colorblind() +
  ggplot2::xlab("Factor") +
  ggplot2::ylab("Average"~eta~"across items")
ggplot2::ggsave(
  plot = plots$exploratory[["eta-averaged-6"]],
  filename = "eta-averaged-6.eps",
  path = here::here("figures", "exploratory"),
  device = cairo_ps,
  width = 9.5, height = 5)

### Model 8 ----
plots$exploratory[["eta-averaged-8"]] <- pars |>
  dplyr::filter(startsWith(par, "eta"), model == 8) |>
  tidyr::separate(par, into = c("par", "item", "class"), extra = "drop") |>
  dplyr::mutate(
    item = factor(as.integer(item), labels = sprintf("par_%s", unique(as.integer(item)))),
    class = as.factor(class)
  ) |>
  dplyr::left_join(item_dictionary, by = "item") |>
  dplyr::mutate(
    item = as.integer(item) |> factor(labels = 1:24)
  ) |>
  dplyr::group_by(factor, class) |>
  dplyr::summarise(mean = mean(estimate)) |>
  ggplot2::ggplot(ggplot2::aes(x = factor, y = mean, group = class, col = class)) +
  ggplot2::geom_line(linewidth = 1) +
  ggthemes::scale_color_colorblind() +
  ggplot2::xlab("Factor") +
  ggplot2::ylab("Average"~eta~"across items")
ggplot2::ggsave(
  plot = plots$exploratory[["eta-averaged-8"]],
  filename = "eta-averaged-8.eps",
  path = here::here("figures", "exploratory"),
  device = cairo_ps,
  width = 9.5, height = 5)


fits[[6]]$par$class_membership |> apply(1, which.max) |> table(classification)
fits[[8]]$par$class_membership |> apply(1, which.max) |> table(classification)
