Sys.setlocale("LC_CTYPE", "C.UTF-8")
# ============================================================================
# EC7062 Business Analytics: figures for the Data 4 (CPS 1991) report
# Tool: R (ggplot2). Run from the R folder:  Rscript EC7062_figures.R
# Inputs: ../data 4.dta (as supplied) and inputs/h5_mret.csv, inputs/coef_cmp.csv
#         (maximum likelihood Heckman results exported from EC7062_analysis.do)
# Output: ../Figures/fig2, fig3, fig5, fig6 and fig7 (PNG, 300 dpi)
# ============================================================================
suppressPackageStartupMessages({
  library(haven); library(ggplot2); library(patchwork); library(quantreg)
  library(sandwich); library(lmtest); library(scales)
})
data_path <- if (file.exists("../data 4.dta")) "../data 4.dta" else "../../data 4.dta"
d <- as.data.frame(zap_labels(read_dta(data_path)))
d$educsq <- d$educ^2
w <- subset(d, inlf == 1)
dir.create("../Figures", showWarnings = FALSE)
dir.create("../Chart_Data", showWarnings = FALSE)
out_csv <- function(x, f) write.csv(x, file.path("../Chart_Data", f), row.names = FALSE)

TEAL <- "#1F4E79"; CORAL <- "#C55A11"; SAND <- "#E9C46A"; INK <- "#000000"; GREY <- "#7F7F7F"
theme_report <- function() {
  theme_minimal(base_size = 10, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", size = 10.5, colour = INK),
          plot.subtitle = element_text(size = 8.5, colour = "#000000"),
          axis.title = element_text(size = 9, colour = INK),
          axis.text = element_text(size = 8.5, colour = "#000000"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.35),
          axis.line = element_line(colour = "#595959", linewidth = 0.35),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8.5, colour = "#000000"),
          plot.margin = margin(6, 10, 6, 6))
}
save_fig <- function(p, file, w, h) ggsave(file.path("../Figures", file), p, width = w, height = h, dpi = 300, bg = "white", type = "cairo")

# ---------------------------------------------------------------- Figure 2
pr <- glm(inlf ~ educ + exper + I(exper^2) + black + hispanic + kidlt6 + kidge6 + I(husearns / 100),
          family = binomial(link = "probit"), data = d)
grid <- expand.grid(educ = 6:18, kidlt6 = c(0, 1))
m <- colMeans(d[, c("exper", "black", "hispanic", "kidge6", "husearns")])
grid$exper <- m["exper"]; grid$black <- m["black"]; grid$hispanic <- m["hispanic"]
grid$husearns <- m["husearns"]; grid$kidge6 <- ifelse(grid$kidlt6 == 1, 0, m["kidge6"])
pp <- predict(pr, newdata = grid, type = "link", se.fit = TRUE)
grid$p <- pnorm(pp$fit); grid$lo <- pnorm(pp$fit - 1.96 * pp$se.fit); grid$hi <- pnorm(pp$fit + 1.96 * pp$se.fit)
grid$grp <- factor(ifelse(grid$kidlt6 == 1, "Child under six", "No child under six"), levels = c("No child under six", "Child under six"))
obs <- aggregate(inlf ~ educ + kidlt6, data = subset(d, educ >= 6), FUN = mean)
obs$n <- aggregate(inlf ~ educ + kidlt6, data = subset(d, educ >= 6), FUN = length)$inlf
obs$grp <- factor(ifelse(obs$kidlt6 == 1, "Child under six", "No child under six"), levels = levels(grid$grp))
p3a <- ggplot(grid, aes(educ, p, colour = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(data = subset(obs, n >= 20), aes(educ, inlf, size = n), shape = 21, colour = "white", alpha = 0.85) +
  scale_colour_manual(values = c(TEAL, CORAL)) + scale_fill_manual(values = c(TEAL, CORAL)) +
  scale_size_area(max_size = 5, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(6, 18, 2)) +
  labs(title = "(a) Probability of participation", subtitle = "Probit prediction at sample means (line, 95% band); observed rates (bubbles, n ≥ 20)",
       x = "Years of schooling", y = "Share in the labour force") + theme_report()
e7 <- pmax(w$educ, 7)
agg <- do.call(rbind, lapply(sort(unique(e7)), function(e) { y <- w$lwage[e7 == e]; data.frame(educ = e, mu = mean(y), se = sd(y) / sqrt(length(y)), n = length(y)) }))
g <- data.frame(educ = seq(7, 18, 0.1))
f1 <- lm(lwage ~ educ, w); f2 <- lm(lwage ~ educ + educsq, w)
g$lin <- predict(f1, transform(g, educsq = educ^2)); g$quad <- predict(f2, transform(g, educsq = educ^2))
p3b <- ggplot(agg, aes(educ, mu)) +
  geom_line(data = g, aes(educ, lin, linetype = "Linear fit"), colour = GREY, linewidth = 0.8) +
  geom_line(data = g, aes(educ, quad, linetype = "Quadratic fit"), colour = TEAL, linewidth = 1) +
  geom_pointrange(aes(ymin = mu - 1.96 * se, ymax = mu + 1.96 * se), colour = INK, size = 0.3, linewidth = 0.5) +
  scale_linetype_manual(values = c("Linear fit" = "dashed", "Quadratic fit" = "solid")) +
  scale_x_continuous(breaks = seq(7, 18, 1), labels = c("≤7", 8:18)) +
  labs(title = "(b) Mean log hourly wage, employed women", subtitle = "Group means with 95% confidence intervals",
       x = "Years of schooling", y = "Log hourly wage") + theme_report()
save_fig(p3a + p3b, "fig2_participation_wages.png", 10, 4.1)
out_csv(grid[, c("educ", "kidlt6", "p", "lo", "hi")], "fig2a_probit_predictions.csv")
out_csv(obs, "fig2a_observed_participation.csv")
out_csv(agg, "fig2b_mean_log_wage_by_schooling.csv")

# ---------------------------------------------------------------- Figure 3
w$grp <- cut(w$educ, c(-1, 11, 12, 15, 18), labels = c("Fewer than 12", "12", "13 to 15", "16 or more"))
lab <- aggregate(lwage ~ grp, data = w, FUN = function(x) c(n = length(x), sd = sd(x)))
lab <- data.frame(grp = lab$grp, n = lab$lwage[, "n"], sd = lab$lwage[, "sd"])
p4 <- ggplot(w, aes(grp, lwage)) +
  geom_violin(fill = TEAL, alpha = 0.18, colour = TEAL, linewidth = 0.4, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.16, outlier.size = 0.6, outlier.colour = CORAL, fill = "white", colour = INK, linewidth = 0.4) +
  geom_text(data = lab, aes(grp, 5.6, label = sprintf("n = %s\nSD = %.2f", comma(n), sd)), size = 2.8, colour = INK, lineheight = 0.9) +
  coord_cartesian(ylim = c(-0.5, 6)) +
  labs(title = "Log hourly wage by schooling group", subtitle = "Violin shows the density; box shows the median and interquartile range; points are outside 1.5 IQR",
       x = "Years of schooling", y = "Log hourly wage") + theme_report()
save_fig(p4, "fig3_wage_distribution.png", 7, 4.3)
q <- do.call(rbind, lapply(levels(w$grp), function(g) { y <- w$lwage[w$grp == g]; data.frame(group = g, n = length(y), mean = mean(y), sd = sd(y), p25 = quantile(y, 0.25), median = median(y), p75 = quantile(y, 0.75)) }))
out_csv(q, "fig3_log_wage_by_group_summary.csv")

# ---------------------------------------------------------------- Figure 5
cc <- read.csv("inputs/coef_cmp.csv")
cc$model <- factor(cc$model, levels = rev(cc$model))
cc$type <- ifelse(grepl("Heckman", cc$model), "Heckman selection model", "OLS, employed women")
p5a <- ggplot(cc, aes(b * 100, model, colour = type)) +
  geom_vline(xintercept = cc$b[cc$model == "OLS M3"] * 100, linetype = "dotted", colour = GREY) +
  geom_errorbar(aes(xmin = lo * 100, xmax = hi * 100), width = 0.18, linewidth = 0.7, orientation = "y") +
  geom_point(size = 2.6) + scale_colour_manual(values = c(CORAL, TEAL)) +
  labs(title = "(a) Linear return across estimators", subtitle = "Point estimate and 95% confidence interval",
       x = "Log points per year of schooling", y = NULL) + theme_report()
m5 <- lm(lwage ~ educ + educsq + exper + I(exper^2) + black + hispanic, w)
V <- vcovHC(m5, type = "HC1"); b <- coef(m5)
mr <- data.frame(e = 8:18)
mr$b <- b["educ"] + 2 * mr$e * b["educsq"]
mr$se <- sqrt(V["educ", "educ"] + 4 * mr$e^2 * V["educsq", "educsq"] + 4 * mr$e * V["educ", "educsq"])
hm <- read.csv("inputs/h5_mret.csv")
p5b <- ggplot(mr, aes(e, b * 100)) +
  geom_ribbon(aes(ymin = (b - 1.96 * se) * 100, ymax = (b + 1.96 * se) * 100), fill = TEAL, alpha = 0.15) +
  geom_line(aes(colour = "OLS M5 (95% band)"), linewidth = 1) +
  geom_line(data = hm, aes(e, b * 100, colour = "Heckman ML, same outcome equation"), linetype = "dashed", linewidth = 0.9) +
  geom_point(data = hm, aes(e, b * 100, colour = "Heckman ML, same outcome equation"), size = 1.8) +
  scale_colour_manual(values = c("OLS M5 (95% band)" = TEAL, "Heckman ML, same outcome equation" = CORAL)) +
  scale_x_continuous(breaks = seq(8, 18, 2)) +
  labs(title = "(b) Marginal return by level of schooling", subtitle = "β₁ + 2β₂s from the quadratic specification",
       x = "Years of schooling (s)", y = "Log points per extra year") + theme_report()
save_fig(p5a + p5b + plot_layout(widths = c(1, 1.1)), "fig5_returns.png", 10, 4.2)
out_csv(cc[, c("model", "b", "lo", "hi")], "fig5a_linear_returns.csv")
out_csv(merge(mr, setNames(hm[, c("e", "b", "se")], c("e", "b_heckman_ml", "se_heckman_ml")), all.x = TRUE), "fig5b_marginal_returns.csv")

# ---------------------------------------------------------------- Figure 6
dd <- data.frame(fit = fitted(m5), res = resid(m5))
dd$z <- (dd$res - mean(dd$res)) / sd(dd$res)
p6a <- ggplot(dd, aes(fit, res)) + geom_hline(yintercept = 0, colour = GREY) +
  geom_point(alpha = 0.25, size = 0.7, colour = TEAL) +
  geom_smooth(method = "loess", se = FALSE, colour = CORAL, linewidth = 0.9, formula = y ~ x) +
  labs(title = "(a) Residuals against fitted values", subtitle = "Loess smoother in orange", x = "Fitted log wage", y = "Residual") + theme_report()
p6b <- ggplot(dd, aes(sample = z)) + stat_qq(size = 0.7, colour = TEAL, alpha = 0.5) + stat_qq_line(colour = CORAL, linewidth = 0.9) +
  labs(title = "(b) Normal quantile plot", subtitle = "Standardised residuals", x = "Theoretical quantile", y = "Sample quantile") + theme_report()
cook <- cooks.distance(m5)
p6c <- ggplot(data.frame(i = seq_along(cook), c = cook), aes(i, c)) +
  geom_segment(aes(xend = i, yend = 0), colour = TEAL, linewidth = 0.3) +
  geom_hline(yintercept = 4 / length(cook), linetype = "dashed", colour = CORAL) +
  labs(title = "(c) Cook's distance", subtitle = "Dashed line: 4/n", x = "Observation", y = "Cook's D") + theme_report()
save_fig(p6a + p6b + p6c, "fig6_diagnostics.png", 10.5, 3.7)
out_csv(data.frame(fitted = dd$fit, residual = dd$res, std_residual = dd$z, cooks_d = cooks.distance(m5)), "fig6_m5_residuals.csv")
cat(sprintf("M5 max Cook's D = %.3f\n", max(cook)))

# ---------------------------------------------------------------- Figure 7
qs <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
qr <- do.call(rbind, lapply(qs, function(q) {
  f <- rq(lwage ~ educ + exper + I(exper^2) + black + hispanic, tau = q, data = w)
  s <- summary(f, se = "nid")$coefficients["educ", ]
  data.frame(q = q, b = s[1], se = s[2])
}))
m3 <- lm(lwage ~ educ + exper + I(exper^2) + black + hispanic, w)
ob <- coef(m3)["educ"]; ose <- sqrt(vcovHC(m3, type = "HC1")["educ", "educ"])
p7 <- ggplot(qr, aes(q, b * 100)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = (ob - 1.96 * ose) * 100, ymax = (ob + 1.96 * ose) * 100, fill = CORAL, alpha = 0.12) +
  geom_hline(yintercept = ob * 100, colour = CORAL, linetype = "dashed") +
  geom_ribbon(aes(ymin = (b - 1.96 * se) * 100, ymax = (b + 1.96 * se) * 100), fill = TEAL, alpha = 0.18) +
  geom_line(colour = TEAL, linewidth = 1) + geom_point(colour = TEAL, size = 2) +
  annotate("text", x = 0.93, y = ob * 100 - 0.55, label = "OLS (M3)", colour = CORAL, size = 3, hjust = 1) +
  scale_x_continuous(breaks = qs, labels = paste0(qs * 100, "th")) +
  labs(title = "Return to schooling across the conditional wage distribution",
       subtitle = "Quantile regression with M3 controls, 95% band; OLS mean effect shaded",
       x = "Quantile of log hourly wage", y = "Log points per year of schooling") + theme_report()
save_fig(p7, "fig7_quantile_process.png", 7, 4.2)
out_csv(qr, "fig7_quantile_estimates.csv")
print(qr)
