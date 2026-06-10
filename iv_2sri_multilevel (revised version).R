################################################################################
# IV and 2SRI for multilevel models 
  
# R code translated based on original Stata code by Eliot Weinstein eliotw@uchicago.edu
# Examples based on my own work, where I'm studying "authors" and "physicians"
# And I have outcome variables that are highly skewed but not suitable for zero-inflated poisson, hence my use of negative binomial
# (with my primary outcome being inpatient length of stay, "LOSCalc"

# Like many students who originally learned and practiced multilevel models in biostats courses
# or courses cross-listed between biostat and regular stat, I learned how to do this in Stata first,
# and then learned how to use for for these tasks later.

# Stata -> R mapping of the main commands used below:
#   regress            -> lm()
#   ivregress 2sls     -> AER::ivreg()
#   poisson            -> glm(family = poisson)
#   logit              -> glm(family = binomial)
#   melogit            -> lme4::glmer(family = binomial)        (multilevel logit)
#   menbreg            -> glmmTMB::glmmTMB(family = nbinom2)    (multilevel neg. binomial)
#   predict, xb        -> predict(fit)                          (linear predictor / fitted)
#   predict, mu / pr   -> predict(fit, type = "response")       (predicted probability/mean)
#   predict, resid     -> residuals(fit)
#   , robust           -> sandwich::vcovHC(type = "HC1") with lmtest::coeftest()
#   vce(cluster id)    -> sandwich::vcovCL(cluster = ~id)
#   irr                -> exp(coef) = incidence-rate ratios
#   estat icc          -> performance::icc() (or compute by hand from VarCorr)
#   coefplot           -> broom(.mixed)::tidy() + ggplot2
#
# Stata factor-variable notation:
#   i.var      -> factor(var)                         (base = lowest level)
#   ib12.var   -> relevel(factor(var), ref = "12")    ("b#" sets the base level)
#
# 2SRI (two-stage residual inclusion): fit a first-stage model for the
# endogenous regressor, take its residual, and include BOTH the endogenous
# regressor and that residual as covariates in the second-stage model. A
# significant residual coefficient is evidence of endogeneity.
#
# Reference for 2SRI in nonlinear models:
#   Terza JV, Basu A, Rathouz PJ (2008). "Two-stage residual inclusion
#   estimation: addressing endogeneity in health econometric modeling."
#   Journal of Health Economics 27(3):531-543.
# https://www.sciencedirect.com/science/article/pii/S0167629607001063
#
# Key point from that paper: for a nonlinear second stage, 2SRI (include the
# residual) is consistent under conditions where two-stage *predictor*
# substitution (plug in the predicted value instead of the observed regressor)
# is NOT. So keep ICU itself in the second stage and add the residual; do not
# replace ICU with its fitted value.
#
# NOTE on standard errors: including a first-stage residual makes the naive
# second-stage SEs too small (they treat the residual as known), and these data
# are clustered. Use the cluster bootstrap in the APPENDIX for inference, not the
# default SEs printed by summary().
################################################################################


# ------------------------------------------------------------------ #
# 0. Packages
# ------------------------------------------------------------------ #
# install.packages(c("haven","lme4","glmmTMB","AER","sandwich","lmtest",
#                     "broom","broom.mixed","ggplot2","dplyr","performance"))

library(haven)        # read .dta
library(dplyr)        # data wrangling
library(lme4)         # glmer (multilevel logit)
library(glmmTMB)      # multilevel negative binomial (menbreg equivalent)
library(AER)          # ivreg (2SLS)
library(sandwich)     # robust / clustered vcov
library(lmtest)       # coeftest
library(broom)        # tidy lm / glm / ivreg
library(broom.mixed)  # tidy mixed models
library(ggplot2)      # plotting
# library(performance)  # icc(); optional helper

# Small helper to build a formula from string pieces (used throughout).
mkform <- function(lhs, rhs) {
  rhs <- rhs[nzchar(rhs)]
  as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))
}


################################################################################
# PART 1 -- GENERIC IV / 2SRI TEMPLATE  (placeholder variable names)
# In the Stata file these used locals: `endog', `instr', `outcome', `controls'.
################################################################################

endog    <- "mean_LOS_MJ1"   # endogenous regressor
instr    <- "instr_var"      # instrument
outcome  <- "outcome_var"    # outcome
controls <- c("cov1", "cov2")# covariates (use character(0) for none)

dat <- read_dta("predicted_LOS_by_author_menbreg.dta")


## 1A. First stage: scatter + linear regression --------------------------------
fs_simple <- lm(mkform(endog, instr), data = dat)
summary(fs_simple)
# For a single instrument the "partial F" is just t^2 on the instrument; a large
# F => strong instrument. (This linear F heuristic applies to a LINEAR first
# stage only -- see the note in Part 2 for the multilevel-logit first stage.)

ggplot(dat, aes(x = .data[[instr]], y = .data[[endog]])) +
  geom_point(color = "navy", size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, color = "maroon", linewidth = 1) +
  labs(title = paste0("First stage: ", endog, " vs ", instr),
       x = paste0("Instrument (", instr, ")"),
       y = paste0("Endogenous regressor (", endog, ")")) +
  theme_minimal()


## 1B. First-stage residual diagnostics ----------------------------------------
dat$endog_hat   <- as.numeric(predict(fs_simple))
dat$endog_resid <- as.numeric(residuals(fs_simple))

ggplot(dat, aes(x = endog_resid)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5, fill = "grey70", color = "white") +
  stat_function(fun = dnorm,
                args = list(mean = mean(dat$endog_resid, na.rm = TRUE),
                            sd   = sd(dat$endog_resid,   na.rm = TRUE)),
                color = "maroon", linewidth = 1) +
  labs(title = "Histogram of first-stage residuals", x = "Residuals", y = "Density") +
  theme_minimal()

ggplot(dat, aes(x = endog_hat, y = endog_resid)) +
  geom_point(size = 1, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "grey60") +
  labs(title = "First-stage residuals vs predicted",
       x = paste0("Predicted ", endog), y = "Residuals") +
  theme_minimal()


## 1C. 2SLS vs OLS (linear outcome) --------------------------------------------
iv_form <- as.formula(paste(outcome, "~",
                            paste(c(endog, controls), collapse = " + "), "|",
                            paste(c(instr, controls), collapse = " + ")))
iv_model  <- ivreg(iv_form, data = dat)
ols_model <- lm(mkform(outcome, c(endog, controls)), data = dat)

coeftest(iv_model,  vcov = vcovHC(iv_model,  type = "HC1"))   # ", robust"
coeftest(ols_model, vcov = vcovHC(ols_model, type = "HC1"))

cmp <- bind_rows(
  tidy(ols_model, conf.int = TRUE) %>% mutate(model = "OLS"),
  tidy(iv_model,  conf.int = TRUE) %>% mutate(model = "IV (2SLS)")
) %>% filter(term == endog)

ggplot(cmp, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  labs(title = paste0("OLS vs IV: ", endog), x = "Coefficient", y = NULL) +
  theme_minimal()


## 1D. 2SRI for a nonlinear (Poisson) outcome ----------------------------------
fs_lin <- lm(mkform(endog, c(instr, controls)), data = dat)
dat$endog_resid <- as.numeric(residuals(fs_lin))

pois_2sri <- glm(mkform(outcome, c(endog, "endog_resid", controls)),
                 data = dat, family = poisson())
coeftest(pois_2sri, vcov = vcovHC(pois_2sri, type = "HC1"))
# Coefficient on endog_resid tests for endogeneity (significant => endogeneity).
# Binary outcome alternative:
# glm(mkform(outcome, c(endog, "endog_resid", controls)), data = dat,
#     family = binomial())


################################################################################
# PART 2 -- MAIN MULTILEVEL 2SRI ANALYSIS
#   First stage : multilevel logit for ICU      (melogit -> glmer)
#   Second stage: multilevel neg. binomial LOS   (menbreg -> glmmTMB), 2SRI, IRR
#   Instrument  : MYHighJeo ;  random intercept by AuthorCode
#
# This part also (a) builds the first-stage residual BOTH ways -- conditional on
# the random effect and at the population level -- and (b) fits the second stage
# as both negative binomial and Poisson, so you can see how the ICU estimate
# moves across those choices. For these data NB is the right family on subject-
# matter grounds (the outcome has no zeros, no fractions, no negatives, and an
# unusual skew/dispersion); Poisson is shown only for comparison.
################################################################################

full <- read_dta("path/to/your_full_dataset.dta")   # <-- edit path

## 2A. Sanity checks ------------------------------------------------------------
req_vars <- c("LOSCalc", "ICU", "MYHighJeo", "AuthorCode")
missing  <- setdiff(req_vars, names(full))
if (length(missing)) warning("Missing variable(s): ", paste(missing, collapse = ", "))

## 2B. Factors with the same base levels as Stata's "ib#." notation ------------
#   ib12.DisD, ib12.Service -> "12";  ib4.MDC -> "4";  ib3.FClass_c -> "3";
#   i.zipCat / i.DayofWeek -> default base (lowest level).
full <- full %>%
  mutate(
    DisD       = relevel(factor(DisD),     ref = "12"),
    Service    = relevel(factor(Service),  ref = "12"),
    MDC        = relevel(factor(MDC),      ref = "4"),
    FClass_c   = relevel(factor(FClass_c), ref = "3"),
    zipCat     = factor(zipCat),
    DayofWeek  = factor(DayofWeek),
    AuthorCode = factor(AuthorCode)
  )

# NOTE on Age: here Age is a nuisance control, so its scale does not matter. If
# in your analysis age is theoretically/causally of interest -- or if you hit
# convergence warnings -- rescale it (e.g. full$Age_z <- scale(full$Age), or use
# age per decade) for interpretable coefficients and better-behaved optimization,
# and substitute the rescaled version into the model formulas below.

controls_ml <- c("DisD", "Age", "sex", "zipCat", "Service", "DayofWeek",
                 "MDC", "FClass_c", "HMProviderCount", "HMServiceCount",
                 "NonHM1", "NonICUNonHM")

## 2C. First stage: multilevel logit (melogit equivalent) ----------------------
fs_form    <- mkform("ICU", c("MYHighJeo", controls_ml, "(1 | AuthorCode)"))
fs_melogit <- glmer(fs_form, data = full, family = binomial,
                    control = glmerControl(optimizer = "bobyqa"))
summary(fs_melogit)

# RELEVANCE / WEAK-INSTRUMENT NOTE (for readers without a causal-inference
# background): with a LINEAR first stage you would check the instrument's partial
# F (rule of thumb F > 10). That rule does NOT transfer to this NONLINEAR
# (logit) first stage. Here, judge instrument relevance from the coefficient and
# z-statistic on MYHighJeo in the summary above (it should be clearly nonzero),
# together with a subject-matter argument that MYHighJeo affects ICU but has no
# direct effect on LOS except through ICU (the exclusion restriction -- which no
# statistical test can verify for a just-identified model like this one).

## 2D. First-stage residual, BOTH ways -----------------------------------------
# Stata's `predict, mu` after melogit returns the prediction INCLUDING the
# empirical-Bayes random effect (conditional). That is re.form = NULL in glmer.
# The population-level prediction (re.form = NA) excludes the random effect; the
# argument for using it is that you do not want the first-stage residual to soak
# up the same author-level heterogeneity that the second-stage random intercept
# is already modeling. Both are defensible -- compute each and check the ICU
# estimate is stable across them.
full$pICU_cond <- as.numeric(predict(fs_melogit, type = "response", re.form = NULL))
full$pICU_marg <- as.numeric(predict(fs_melogit, type = "response", re.form = NA))
full$rICU_cond <- full$ICU - full$pICU_cond   # conditional residual (matches Stata)
full$rICU_marg <- full$ICU - full$pICU_marg   # population-level residual

## 2E. Second stage: NB and Poisson, conditional and marginal residual ---------
ss_terms <- c("ICU", "MYHighJeo", controls_ml, "(1 | AuthorCode)")

# Primary specification: negative binomial, conditional residual.
ss_nb_cond  <- glmmTMB(mkform("LOSCalc", c(ss_terms, "rICU_cond")),
                       data = full, family = nbinom2)
# Same, population-level residual.
ss_nb_marg  <- glmmTMB(mkform("LOSCalc", c(ss_terms, "rICU_marg")),
                       data = full, family = nbinom2)
# Poisson comparison (conditional residual).
ss_poi_cond <- glmmTMB(mkform("LOSCalc", c(ss_terms, "rICU_cond")),
                       data = full, family = poisson)

summary(ss_nb_cond)

# Incidence-rate ratios for the primary model (Stata: , irr).
tidy(ss_nb_cond, effects = "fixed", exponentiate = TRUE, conf.int = TRUE)

# Compare the ICU coefficient across family x residual choices.
ss_compare <- bind_rows(
  tidy(ss_nb_cond,  effects = "fixed") %>% mutate(spec = "NB, conditional resid"),
  tidy(ss_nb_marg,  effects = "fixed") %>% mutate(spec = "NB, marginal resid"),
  tidy(ss_poi_cond, effects = "fixed") %>% mutate(spec = "Poisson, conditional resid")
) %>% filter(term == "ICU") %>% select(spec, estimate, std.error)
print(ss_compare)

## 2F. Diagnostics -------------------------------------------------------------
# Endogeneity test: is the first-stage residual significant in the second stage?
summary(ss_nb_cond)$coefficients$cond["rICU_cond", , drop = FALSE]

# NB vs Poisson on these data: compare AIC (lower = better). NB should win
# decisively given the dispersion; this is just a demonstration of the check.
AIC(ss_nb_cond, ss_poi_cond)

# ICC for the first-stage multilevel logit (Stata: estat icc), latent scale.
# performance::icc(fs_melogit)   # convenient, or by hand:
vc_re   <- as.data.frame(VarCorr(fs_melogit)$AuthorCode)[1, 1]
icc_log <- vc_re / (vc_re + pi^2 / 3)
cat("First-stage (logit) ICC:", round(icc_log, 4), "\n")


################################################################################
# PART 3 -- TOP-20 vs FULL-SAMPLE COMPARISON
# Keep the two samples as separate data frames and refit. Replace the filter
# with your real Top-20 definition.
################################################################################

top20_ids <- full %>% count(AuthorCode, sort = TRUE) %>% slice_head(n = 20) %>%
  pull(AuthorCode)
top20 <- full %>% filter(AuthorCode %in% top20_ids) %>%
  mutate(AuthorCode = droplevels(AuthorCode))

fs_top20 <- glmer(fs_form, data = top20, family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))
top20$rICU_cond <- top20$ICU - as.numeric(predict(fs_top20, type = "response",
                                                  re.form = NULL))
ss_top20 <- glmmTMB(mkform("LOSCalc", c(ss_terms, "rICU_cond")),
                    data = top20, family = nbinom2)

comparison <- bind_rows(
  tidy(ss_nb_cond, effects = "fixed", conf.int = TRUE) %>% mutate(sample = "Full"),
  tidy(ss_top20,   effects = "fixed", conf.int = TRUE) %>% mutate(sample = "Top-20")
) %>% select(sample, term, estimate, std.error, p.value, conf.low, conf.high)
print(comparison, n = Inf)


################################################################################
# PART 4 -- FOREST PLOT: ICU effect, "Regular" (naive) vs 2SRI
################################################################################

ss_naive <- glmmTMB(mkform("LOSCalc", ss_terms), data = full, family = nbinom2)

forest_df <- bind_rows(
  tidy(ss_naive,   effects = "fixed", conf.int = TRUE) %>% mutate(model = "Regular"),
  tidy(ss_nb_cond, effects = "fixed", conf.int = TRUE) %>% mutate(model = "2SRI")
) %>% filter(term == "ICU")

p_forest <- ggplot(forest_df, aes(x = estimate, y = model, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_color_manual(values = c("Regular" = "blue", "2SRI" = "red")) +
  labs(title = "ICU treatment effect: Regular vs 2SRI",
       x = "Coefficient (log scale)", y = NULL, color = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
print(p_forest)

ggsave("forest_pub.png",    p_forest, width = 12, height = 8, units = "in", dpi = 100)
ggsave("forest_poster.png", p_forest, width = 36, height = 48, units = "in", dpi = 150,
       limitsize = FALSE)

write.csv(
  forest_df %>% transmute(model, estimate, std.error, conf.low, conf.high, p.value),
  "forest_table.csv", row.names = FALSE)


################################################################################
# PART 5 -- ALL-COEFFICIENT FOREST PLOT  (LPM-first-stage workflow)
#   naive multilevel NB; LPM first stage for ICU on Z1 Z2; 2SRI second stage;
#   extract every non-intercept coefficient from both models and plot.
################################################################################

d <- read_dta("path/to/your_full_dataset.dta")   # <-- edit path
d$physician_id <- factor(d$physician_id)
covs2 <- c("age", "male", "comorbidity_score", "othercov1", "othercov2")

m_naive <- glmmTMB(mkform("LOS", c("ICU", covs2, "(1 | physician_id)")),
                   data = d, family = nbinom2)

fs_lpm <- lm(mkform("ICU", c("Z1", "Z2", covs2)), data = d)
coeftest(fs_lpm, vcov = vcovCL(fs_lpm, cluster = ~ physician_id))
d$ICU_resid <- as.numeric(residuals(fs_lpm))
# Probit first-stage alternative:
# fs_probit <- glm(mkform("ICU", c("Z1","Z2",covs2)), data = d,
#                  family = binomial(link = "probit"))
# d$ICU_resid <- d$ICU - predict(fs_probit, type = "response")

m_2sri <- glmmTMB(mkform("LOS", c("ICU", "ICU_resid", covs2, "(1 | physician_id)")),
                  data = d, family = nbinom2)

N_obs      <- nobs(m_2sri)
N_clusters <- nlevels(droplevels(model.frame(m_2sri)$physician_id))
cat("Observations:", N_obs, " | Physician clusters:", N_clusters, "\n")

extract_coefs <- function(model, label) {
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    transmute(model = label, term, estimate,
              lower = conf.low, upper = conf.high, p = p.value)
}
all_coefs <- bind_rows(
  extract_coefs(m_naive, "Naive_menbreg"),
  extract_coefs(m_2sri,  "2SRI_menbreg")
) %>%
  mutate(est_ci = sprintf("%.3f (%.3f, %.3f)", estimate, lower, upper)) %>%
  arrange(term, model)

datestr <- format(Sys.Date(), "%Y-%m-%d")
write.csv(all_coefs, paste0("forest_table_recreated_", datestr, ".csv"),
          row.names = FALSE)
write.csv(data.frame(N = N_obs, clusters = N_clusters),
          paste0("forest_summary_recreated_", datestr, ".csv"), row.names = FALSE)
print(all_coefs %>% select(model, term, est_ci, p), n = Inf)

# If the factor controls expand to many dummies this plot gets unreadable;
# filter `all_coefs` to the terms you actually want to show before plotting.
p_all <- ggplot(all_coefs, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("Naive_menbreg" = "blue", "2SRI_menbreg" = "red")) +
  labs(title = "Forest plot: naive vs 2SRI (all coefficients)",
       x = "Coefficient", y = NULL, color = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
print(p_all)

ggsave(paste0("forest_pub_recreated_",    datestr, ".png"),
       p_all, width = 12, height = 8, dpi = 100)
ggsave(paste0("forest_poster_recreated_", datestr, ".png"),
       p_all, width = 36, height = 48, dpi = 100, limitsize = FALSE)


################################################################################
# APPENDIX -- CLUSTER BOOTSTRAP for 2SRI standard errors  (use this for inference)
#
# Why: the second-stage SEs printed above are too small because they treat the
# first-stage residual as known, and the data are clustered. The fix is a
# nonparametric CLUSTER bootstrap: resample whole clusters with replacement,
# refit BOTH stages inside each replicate, and collect the ICU coefficient.
#
# Two important details:
#  (1) Resample at the TOP level of clustering. If authors/physicians are the top
#      level, resample AuthorCode (cluster_boot_flat). If authors are nested
#      within departments and the dependence lives at the department level,
#      resample dept and carry the authors along (cluster_boot_nested).
#  (2) When a cluster is drawn more than once, its copies must be treated as
#      DISTINCT clusters (fresh ids), otherwise the random-effect variance is
#      underestimated. Both functions below relabel copies accordingly.
#
# These refit glmer + glmmTMB hundreds of times and are slow. Start with a small
# R to check it runs, then scale up; parallelize with furrr::future_map if needed.
################################################################################

## A1. Flat case: resample AuthorCode (matches Part 2's (1 | AuthorCode)) -------
cluster_boot_flat <- function(data, cluster_var, fixed_first, fixed_second,
                              y1, y2, family2 = nbinom2, re_pred = NULL,
                              R = 200, seed = 1) {
  set.seed(seed)
  ids  <- unique(data[[cluster_var]])
  ests <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    drawn <- sample(ids, length(ids), replace = TRUE)
    bd <- do.call(rbind, Map(function(id, b) {
      ch <- data[data[[cluster_var]] == id, , drop = FALSE]
      ch$.grp <- b                       # fresh cluster id for this draw
      ch
    }, drawn, seq_along(drawn)))
    bd$.grp <- factor(bd$.grp)
    ests[r] <- tryCatch({
      f1 <- glmer(mkform(y1, c(fixed_first, "(1 | .grp)")),
                  data = bd, family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))
      bd$.resid <- bd[[y1]] - as.numeric(predict(f1, type = "response",
                                                 re.form = re_pred))
      f2 <- glmmTMB(mkform(y2, c(fixed_second, ".resid", "(1 | .grp)")),
                    data = bd, family = family2)
      fixef(f2)$cond[[y1]]
    }, error = function(e) NA_real_)
  }
  ests
}

# Example call (point estimate stays ss_nb_cond from Part 2):
# boot_flat <- cluster_boot_flat(
#   data        = full,
#   cluster_var = "AuthorCode",
#   fixed_first = c("MYHighJeo", controls_ml),     # first-stage RHS (incl instrument)
#   fixed_second= c("ICU", "MYHighJeo", controls_ml), # second-stage RHS (incl ICU)
#   y1 = "ICU", y2 = "LOSCalc", family2 = nbinom2,
#   re_pred = NULL,   # NULL = conditional residual; NA = population-level residual
#   R = 200, seed = 1)
# point <- fixef(ss_nb_cond)$cond[["ICU"]]
# cat("ICU coef:", round(point, 4),
#     " 95% cluster-bootstrap CI:",
#     paste(round(quantile(boot_flat, c(.025, .975), na.rm = TRUE), 4), collapse = " to "),
#     "\n")
# # On the IRR scale: exp(point) and exp(quantile(...)).

## A2. Nested case: authors within departments --------------------------------
# Use when the model is (1 | dept/AuthorCode), i.e. (1 | dept) + (1 | dept:AuthorCode),
# and the top level of dependence is the department. We resample DEPARTMENTS with
# replacement; each drawn department copy gets a fresh dept id, and its authors
# are made unique within that copy so duplicated departments are independent.
cluster_boot_nested <- function(data, dept_var, author_var,
                                fixed_first, fixed_second, y1, y2,
                                family2 = nbinom2, re_pred = NULL,
                                R = 200, seed = 1) {
  set.seed(seed)
  depts <- unique(data[[dept_var]])
  ests  <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    drawn <- sample(depts, length(depts), replace = TRUE)
    bd <- do.call(rbind, Map(function(dp, b) {
      ch <- data[data[[dept_var]] == dp, , drop = FALSE]
      ch$.dept <- b                                   # fresh department id
      ch$.auth <- paste0(b, "_", ch[[author_var]])    # author unique within copy
      ch
    }, drawn, seq_along(drawn)))
    bd$.dept <- factor(bd$.dept)
    bd$.auth <- factor(bd$.auth)
    ests[r] <- tryCatch({
      f1 <- glmer(mkform(y1, c(fixed_first, "(1 | .dept)", "(1 | .auth)")),
                  data = bd, family = binomial,
                  control = glmerControl(optimizer = "bobyqa"))
      bd$.resid <- bd[[y1]] - as.numeric(predict(f1, type = "response",
                                                 re.form = re_pred))
      f2 <- glmmTMB(mkform(y2, c(fixed_second, ".resid",
                                 "(1 | .dept)", "(1 | .auth)")),
                    data = bd, family = family2)
      fixef(f2)$cond[[y1]]
    }, error = function(e) NA_real_)
  }
  ests
}
# Example call mirrors cluster_boot_flat but passes dept_var = "dept",
# author_var = "AuthorCode", and your point-estimate model must likewise use
# (1 | dept) + (1 | dept:AuthorCode).

cat("Done. Tables and figures written to the working directory:", getwd(), "\n")
################################################################################
