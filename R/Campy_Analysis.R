# ------------------------------------------------------------------
# Campylobacteriosis in Switzerland, 2013-2026:
# seasonality, the festive winter peak and the effect of temperature
#
# Data (see README for sources and licences)
#   data/raw/BAG_CAMPYLOBACTERIOSIS_oblig_latest.csv   FOPH/BAG IDD, mandatory reporting
#   data/raw/BAG_SALMONELLOSIS_oblig_latest.csv        comparison disease
#   data/raw/ogd-nbcn_<station>_m.csv                   MeteoSwiss NBCN homogeneous monthly series
#
# Run from the repository root:  Rscript R/Campy_Analysis.R
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(MASS); library(splines)
  library(sandwich); library(ggplot2); library(patchwork); library(scales)
})
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

START <- as.Date("2013-01-01")                 # first month in the BAG export
PANDEMIC <- c(as.Date("2020-03-01"), as.Date("2021-06-01"))   # excluded in a sensitivity analysis
STATIONS <- c(SMA = "Zürich/Fluntern", BAS = "Basel/Binningen", BER = "Bern/Zollikofen",
              GVE = "Genève/Cointrin", LUG = "Lugano", LUZ = "Luzern", STG = "St. Gallen")
HAC_LAG <- 3                                   # Newey-West lag (months) for robust standard errors

# ---------------- 1. Case counts ----------------
read_bag <- function(file, disease) {
  d <- fread(file, na.strings = "NA")
  if (!"type" %in% names(d)) d[, type := "all"]
  m <- d[valueCategory == "cases" & temporal_type == "month" & georegion == "CHFL" &
           agegroup == "all" & sex == "all" & type == "all",
         .(month = as.Date(paste0(sub("-M", "-", temporal), "-01")), cases = as.numeric(value),
           pop = as.numeric(pop), complete = dataComplete)]
  y <- d[valueCategory == "cases" & temporal_type == "year" & georegion == "CHFL" &
           agegroup == "all" & sex == "all" & type == "all",
         .(year = as.integer(temporal), cases = as.numeric(value), pop = as.numeric(pop))]
  stopifnot(!anyDuplicated(m$month), all(!is.na(m$cases)))
  list(month = m[order(month)][, disease := disease], year = y[order(year)][, disease := disease])
}
campy <- read_bag("data/raw/BAG_CAMPYLOBACTERIOSIS_oblig_latest.csv", "Campylobacteriosis")
salmo <- read_bag("data/raw/BAG_SALMONELLOSIS_oblig_latest.csv", "Salmonellosis")
cases <- rbind(campy$month, salmo$month)
annual <- rbind(campy$year, salmo$year)

# consistency check: monthly counts add up to the published annual totals
chk <- merge(cases[, .(sum_months = sum(cases), n_months = .N), by = .(disease, year = year(month))],
             annual[, .(disease, year, annual = cases)], by = c("disease", "year"))
chk[, diff := sum_months - annual]
fwrite(chk, "results/check_monthly_vs_annual.csv")

# ---------------- 2. Temperature ----------------
read_station <- function(code) {
  f <- sprintf("data/raw/ogd-nbcn_%s_m.csv", tolower(code))
  d <- fread(f, sep = ";", na.strings = "", encoding = "Latin-1")
  d[, .(station = code, month = as.Date(substr(reference_timestamp, 1, 10), "%d.%m.%Y"),
        tmean = as.numeric(ths200m0), precip = as.numeric(rhs150m0))]
}
met <- rbindlist(lapply(names(STATIONS), read_station))[month >= START - 31]
met_n <- met[, .(n_t = sum(!is.na(tmean)), n_p = sum(!is.na(precip))), by = month]
stopifnot(all(met_n[month <= max(cases$month)]$n_t == length(STATIONS)))
temp <- met[, .(tmean = mean(tmean), tmean_zh = tmean[station == "SMA"], precip = mean(precip, na.rm = TRUE)), by = month]
setorder(temp, month)
# anomaly = deviation from the 2013-2026 mean of the same calendar month
temp[, cal := month(month)]
temp[month >= START & month <= max(cases$month), `:=`(clim = mean(tmean), clim_zh = mean(tmean_zh), clim_p = mean(precip)), by = cal]
temp[, `:=`(clim = clim[!is.na(clim)][1], clim_zh = clim_zh[!is.na(clim_zh)][1], clim_p = clim_p[!is.na(clim_p)][1]), by = cal]
temp[, `:=`(anom0 = tmean - clim, anom0_zh = tmean_zh - clim_zh, panom0 = precip - clim_p)]
temp[, `:=`(anom1 = shift(anom0), anom1_zh = shift(anom0_zh), tmean1 = shift(tmean))]
fwrite(temp[month >= START - 31], "results/temperature_monthly.csv")

# ---------------- 3. Analysis data ----------------
dat <- merge(cases, temp[, .(month, tmean, tmean1, anom0, anom1, anom0_zh, anom1_zh, panom0)], by = "month")
dat[, `:=`(year = year(month), cal = month(month), days = as.integer(as.Date(format(month + 32, "%Y-%m-01")) - month))]
dat[, `:=`(fyear = factor(year), fcal = factor(cal), jan = as.integer(cal == 1), dec = as.integer(cal == 12),
           pandemic = month >= PANDEMIC[1] & month <= PANDEMIC[2],
           inc_day = cases / pop / days * 1e5)]            # cases per 100,000 per day
dat[, off := log(pop) + log(days)]
stopifnot(!anyNA(dat[, .(cases, pop, tmean, tmean1, anom0, anom1)]))
fwrite(dat[, .(disease, month, cases, pop, days, inc_day, tmean, tmean1, anom0, anom1, pandemic)],
       "results/analysis_data.csv")

# ---------------- helpers ----------------
rr_row <- function(fit, terms, label, V = vcovHAC_nb(fit)) {
  b <- coef(fit); w <- setNames(rep(0, length(b)), names(b)); w[terms] <- 1
  est <- sum(w * b); se <- sqrt(drop(t(w) %*% V %*% w))
  data.table(term = label, rr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
             p = 2 * pnorm(-abs(est / se)))
}
vcovHAC_nb <- function(fit) NeweyWest(fit, lag = HAC_LAG, prewhite = FALSE, adjust = TRUE)
acf1 <- function(fit) { r <- residuals(fit, type = "pearson"); cor(r[-1], r[-length(r)]) }

# ---------------- Q1. Temperature anomalies (primary model) ----------------
# log E[cases] = calendar month + year + b0*anomaly(t) + b1*anomaly(t-1) + log(pop*days)
# Calendar-month and year effects absorb the seasonal cycle, long-term trends and the
# pandemic years; the temperature effect is estimated from months that were warmer or
# colder than usual for the time of year.
f_temp <- cases ~ fcal + fyear + anom0 + anom1 + offset(off)
fit_temp <- function(d, f = f_temp) glm.nb(f, data = d)
q1 <- rbindlist(lapply(split(dat, by = "disease"), function(d) {
  fit <- fit_temp(d)
  r <- rbind(rr_row(fit, "anom0", "same month"), rr_row(fit, "anom1", "previous month"),
             rr_row(fit, c("anom0", "anom1"), "both months (sum)"))
  r[, `:=`(disease = d$disease[1], theta = fit$theta, n = nobs(fit), resid_acf1 = acf1(fit))]
}))
fwrite(q1, "results/q1_temperature_rr.csv")

# sensitivity analyses for campylobacteriosis, combined effect (+1 °C in both months)
dc <- dat[disease == "Campylobacteriosis"]; ds <- dat[disease == "Salmonellosis"]
dc[, log_prev := log(shift(cases) / shift(pop) / shift(days) * 1e5)]
sens <- function(label, fit, terms = c("anom0", "anom1"), V = vcovHAC_nb(fit))
  cbind(analysis = label, rr_row(fit, terms, "both months (sum)", V))
fit_c <- fit_temp(dc)
qp <- glm(cases ~ fcal + fyear + anom0 + anom1 + offset(off), family = quasipoisson, data = dc)
q1s <- rbind(
  sens("Main model (7 stations, HAC SE)", fit_c),
  sens("Model-based SE instead of HAC", fit_c, V = vcov(fit_c)),
  sens("Quasi-Poisson", qp, V = NeweyWest(qp, lag = HAC_LAG, prewhite = FALSE, adjust = TRUE)),
  sens("Pandemic months excluded (Mar 2020-Jun 2021)", fit_temp(dc[pandemic == FALSE])),
  sens("December and January excluded", fit_temp(dc[!cal %in% c(12, 1)])),
  sens("Zürich station only", glm.nb(cases ~ fcal + fyear + anom0_zh + anom1_zh + offset(off), data = dc),
       terms = c("anom0_zh", "anom1_zh")),
  sens("Adjusted for precipitation anomaly", glm.nb(cases ~ fcal + fyear + anom0 + anom1 + panom0 + offset(off), data = dc)),
  sens("Complete years only (2013-2025)", fit_temp(dc[year <= 2025])),
  sens("Smooth time trend instead of year effects", glm.nb(cases ~ fcal + ns(as.numeric(month), df = 28) + anom0 + anom1 + offset(off), data = dc)),
  sens("Adjusted for previous month's cases (autoregressive term)",
       glm.nb(cases ~ fcal + fyear + anom0 + anom1 + log_prev + offset(off), data = dc[!is.na(log_prev)]))
)
fwrite(q1s, "results/q1_temperature_sensitivity.csv")

# ---------------- Q2. Festive winter peak ----------------
# Seasonality is modelled through temperature (natural splines, same and previous month)
# plus year effects. January and December indicators then measure how far these months
# lie above what their temperature and year predict. Salmonellosis serves as a negative
# control: it is also food-borne and reported through the same system, but it is not
# known to be associated with meat fondue.
f_fest <- cases ~ fyear + ns(tmean, df = 3) + ns(tmean1, df = 3) + jan + dec + offset(off)
q2 <- rbindlist(lapply(split(dat, by = "disease"), function(d) {
  fit <- glm.nb(f_fest, data = d)
  r <- rbind(rr_row(fit, "jan", "January"), rr_row(fit, "dec", "December"))
  r[, `:=`(disease = d$disease[1], model = "temperature splines + year", resid_acf1 = acf1(fit))]
  fit2 <- glm.nb(update(f_fest, . ~ . + sin(2 * pi * cal / 12) + cos(2 * pi * cal / 12)), data = d)
  r2 <- rbind(rr_row(fit2, "jan", "January"), rr_row(fit2, "dec", "December"))
  r2[, `:=`(disease = d$disease[1], model = "+ one seasonal harmonic", resid_acf1 = acf1(fit2))]
  fit3 <- glm.nb(f_fest, data = d[pandemic == FALSE])
  r3 <- rbind(rr_row(fit3, "jan", "January"), rr_row(fit3, "dec", "December"))
  r3[, `:=`(disease = d$disease[1], model = "pandemic months excluded", resid_acf1 = acf1(fit3))]
  fit4 <- glm.nb(cases ~ fyear + ns(tmean, df = 5) + ns(tmean1, df = 5) + jan + dec + offset(off), data = d)
  r4 <- rbind(rr_row(fit4, "jan", "January"), rr_row(fit4, "dec", "December"))
  r4[, `:=`(disease = d$disease[1], model = "temperature splines with 5 df", resid_acf1 = acf1(fit4))]
  rbind(r, r2, r3, r4)
}))
fwrite(q2, "results/q2_festive_peak_rr.csv")

# ratio of the January (and December) excess, campylobacteriosis vs salmonellosis
# (log ratio; SEs from the two separate HAC fits combined as if independent)
q2_ratio <- rbindlist(lapply(c("January", "December"), function(tm) {
  a <- q2[model == "temperature splines + year" & term == tm & disease == "Campylobacteriosis"]
  b <- q2[model == "temperature splines + year" & term == tm & disease == "Salmonellosis"]
  se <- function(x) (log(x$hi) - log(x$lo)) / (2 * 1.96)
  est <- log(a$rr) - log(b$rr); s <- sqrt(se(a)^2 + se(b)^2)
  data.table(term = tm, ratio = exp(est), lo = exp(est - 1.96 * s), hi = exp(est + 1.96 * s), p = 2 * pnorm(-abs(est / s)))
}))
fwrite(q2_ratio, "results/q2_campy_vs_salmo_ratio.csv")

# January excess by year and its time trend
jan_year <- rbindlist(lapply(split(dat, by = "disease"), function(d) {
  d <- copy(d); d[, jan_y := factor(ifelse(cal == 1, year, 0))]
  fit <- glm.nb(cases ~ fyear + ns(tmean, df = 3) + ns(tmean1, df = 3) + dec + jan_y + offset(off), data = d)
  V <- vcov(fit)   # one observation per indicator: model-based SE (a HAC estimate would be degenerate)
  rbindlist(lapply(sort(unique(d[cal == 1]$year)), function(y)
    cbind(disease = d$disease[1], year = y, rr_row(fit, paste0("jan_y", y), "January", V))))
}))
# absolute excess: observed January cases minus cases expected without the excess
jan_year <- merge(jan_year, dat[cal == 1, .(disease, year, observed = cases)], by = c("disease", "year"))
jan_year[, `:=`(expected = observed / rr, excess_cases = observed - observed / rr)]
fwrite(jan_year, "results/q2_january_excess_by_year.csv")

jan_trend <- rbindlist(lapply(split(dat, by = "disease"), function(d) {
  d <- copy(d); d[, jan_t := jan * (year - 2019)]
  fit <- glm.nb(cases ~ fyear + ns(tmean, df = 3) + ns(tmean1, df = 3) + dec + jan + jan_t + offset(off), data = d)
  cbind(disease = d$disease[1], rbind(rr_row(fit, "jan", "January excess in 2019 (centre)"),
                                      rr_row(fit, "jan_t", "change in January excess per year")))
}))
fwrite(jan_trend, "results/q2_january_trend.csv")

# ---------------- Q3. Long-term level ----------------
inc_year <- annual[year < 2026, .(disease, year, cases, pop, inc_100k = cases / pop * 1e5)]
fwrite(inc_year, "results/q3_annual_incidence.csv")
ytd <- dat[cal <= max(dat[year == 2026]$cal),
           .(cases = sum(cases), inc_100k = sum(cases) / mean(pop) * 1e5,
             tmean = mean(tmean), anom = mean(anom0)), by = .(disease, year)]
ytd_ref <- ytd[year %in% 2023:2025, .(ref_cases = mean(cases), ref_inc = mean(inc_100k)), by = disease]
ytd <- merge(ytd, ytd_ref, by = "disease")[, ratio_to_2023_25 := inc_100k / ref_inc]
fwrite(ytd, "results/q3_jan_aug_by_year.csv")

# ---------------- Summary numbers used in text ----------------
seas <- dat[year <= 2025, .(mean_cases = mean(cases), mean_inc_day = mean(inc_day)), by = .(disease, cal)]
fwrite(seas, "results/seasonal_profile_2013_2025.csv")
fwrite(data.table(station = names(STATIONS), name = STATIONS), "results/stations.csv")

# input manifest and session info
man <- data.table(file = list.files("data/raw", full.names = TRUE))
man[, `:=`(bytes = file.size(file), sha256 = vapply(file, function(f)
  sub(" .*", "", system2("sha256sum", f, stdout = TRUE)), ""))]
fwrite(man, "results/input_manifest.csv")
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

# ---------------- Figures ----------------
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GREY <- "#c9c7c1"; GREY_D <- "#6f6d68"
TEXT <- "#1f1f1e"; TEXT2 <- "#52514e"; GRID <- "#e8e7e3"
COLS <- c(Campylobacteriosis = BLUE, Salmonellosis = ORANGE)
theme_p <- theme_minimal(base_size = 20) +
  theme(text = element_text(colour = TEXT), axis.text = element_text(colour = TEXT2, size = 17),
        axis.title = element_text(colour = TEXT2, size = 18), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = GRID, linewidth = 0.45),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.text = element_text(size = 17), plot.margin = margin(8, 14, 6, 6))
save <- function(p, name, w = 10, h = 6) ggsave(file.path("figures", name), p, width = w, height = h, dpi = 200, bg = "white")

# A: monthly incidence, January highlighted
pa_dat <- dat[disease == "Campylobacteriosis"]
pA <- ggplot(pa_dat, aes(month, inc_day * 30)) +
  annotate("rect", xmin = PANDEMIC[1], xmax = PANDEMIC[2] + 30, ymin = -Inf, ymax = Inf, fill = "#f1f0ec") +
  annotate("text", x = PANDEMIC[1] + 243, y = 15.6, label = "Pandemic\nmeasures", size = 5, colour = GREY_D, vjust = 1, lineheight = 0.9) +
  geom_line(colour = BLUE, linewidth = 0.9) +
  geom_point(data = pa_dat[cal == 1], colour = TEXT, fill = "white", shape = 21, size = 3.4, stroke = 1.3) +
  annotate("point", x = as.Date("2013-11-01"), y = 15.1, shape = 21, colour = TEXT, fill = "white", size = 3.4, stroke = 1.3) +
  annotate("text", x = as.Date("2014-01-15"), y = 15.1, label = "January", hjust = 0, size = 5.6, colour = TEXT2) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y", expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Cases per 100,000 per 30 days") + theme_p
save(pA, "fig_A_timeseries.png", 10, 5.7)

# B: incidence against temperature; January and December marked
pb_dat <- dat[disease == "Campylobacteriosis"][, grp := fifelse(cal == 1, "January", fifelse(cal == 12, "December", "Other months"))]
pb_dat[, grp := factor(grp, c("Other months", "December", "January"))]
pB <- ggplot(pb_dat, aes(tmean, inc_day * 30)) +
  geom_point(aes(fill = grp, size = grp), shape = 21, colour = "white", stroke = 0.6) +
  scale_fill_manual(values = c("Other months" = "#9fc3ee", December = GREY_D, January = TEXT)) +
  scale_size_manual(values = c("Other months" = 3.3, December = 4, January = 4.4)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Mean temperature of the month (°C, 7 stations)", y = "Cases per 100,000 per 30 days") +
  theme_p + theme(panel.grid.major.x = element_line(colour = GRID, linewidth = 0.45)) +
  guides(size = "none", fill = guide_legend(override.aes = list(size = 5)))
save(pB, "fig_B_temperature_scatter.png", 10, 5.7)

# C: January excess by year, campylobacteriosis vs salmonellosis
pC <- ggplot(jan_year, aes(year + ifelse(disease == "Campylobacteriosis", -0.14, 0.14), rr, colour = disease)) +
  geom_hline(yintercept = 1, colour = GREY_D, linewidth = 0.5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 1.1) +
  geom_point(size = 3.4) +
  scale_colour_manual(values = COLS) +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3), labels = c("0.5", "0.75", "1", "1.5", "2", "3")) +
  scale_x_continuous(breaks = seq(2013, 2026, 2)) +
  labs(x = NULL, y = "January: observed / expected\n(log scale)") + theme_p
save(pC, "fig_C_january_by_year.png", 10, 5.7)

# D: temperature effect
pd <- copy(q1)[, term := factor(term, rev(c("same month", "previous month", "both months (sum)")),
                                  rev(c("Same month", "Previous month", "Both months")))]
pd[, disease := factor(disease, c("Salmonellosis", "Campylobacteriosis"))]
pD <- ggplot(pd, aes(rr, term, colour = disease)) +
  geom_vline(xintercept = 1, colour = GREY_D, linewidth = 0.5) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1.3, position = position_dodge(width = 0.55)) +
  geom_point(size = 3.8, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = COLS, breaks = c("Campylobacteriosis", "Salmonellosis")) +
  scale_x_continuous(labels = function(x) sprintf("%+.0f%%", (x - 1) * 100)) +
  labs(x = "Change in cases per +1 °C warmer than usual", y = NULL) +
  theme_p + theme(panel.grid.major.x = element_line(colour = GRID, linewidth = 0.45), panel.grid.major.y = element_blank())
save(pD, "fig_D_temperature_effect.png", 10, 5.7)

cat("\n=== Q1 temperature ===\n"); print(q1)
cat("\n=== Q1 sensitivity ===\n"); print(q1s)
cat("\n=== Q2 festive peak ===\n"); print(q2)
cat("\n=== Q2 January trend ===\n"); print(jan_trend)
cat("\n=== Q2 ratio ===\n"); print(q2_ratio)
cat("\n=== check monthly vs annual ===\n"); print(chk[diff != 0])
cat("\n=== Jan-Aug ===\n"); print(ytd[year >= 2022])
