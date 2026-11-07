library(quantmod)
library(QRM)
library(tseries)
set.seed(1234)

tickers <- c("MC.PA","OR.PA","SAP.DE","ASML.AS","SIE.DE","TTE.PA","ENI.MI","BNP.PA","ISP.MI","ADS.DE")
sector_map <- data.frame(
  ticker = tickers,
  sector = c("Consumer","Consumer","Technology","Technology","Industrials",
             "Energy","Energy","Financials","Financials","Consumer"))

start_date <- "2017-01-01"
end_date_in_sample <- "2021-12-31"
start_date_out_sample <- "2022-01-01"
end_date <- "2024-12-31"
split.date  <- as.Date(end_date_in_sample)
shock.dates <- as.Date(c("2020-02-24", "2022-02-24"))

price_list <- list()
for (tk in tickers) {
  price_list[[tk]] <- getSymbols(tk, src = "yahoo",
                                 from = start_date, to = end_date,
                                 auto.assign = FALSE)}

adj_list <- lapply(price_list, Ad)
prices <- do.call(merge, adj_list)
colnames(prices) <- names(adj_list)

prices_filled <- na.locf(prices, na.rm = FALSE)
prices_clean <- na.omit(prices_filled)
log_returns <- diff(log(prices_clean))[-1,]

returns_in <- log_returns[paste(start_date, end_date_in_sample, sep = "/")]

par(mfrow = c(2, 1))
plot_series <- function(x, title, ylab, legend = FALSE) {
  d <- index(x); m <- coredata(x)
  matplot(d, m, type = "l", lty = 1, col = 1:ncol(m),
          xlab = "", ylab = ylab, main = title, xaxt = "n")
  axis.Date(1, at = seq(min(d), max(d), by = "year"), format = "%Y")
  abline(v = split.date,  col = "black", lty = 2, lwd = 2)
  abline(v = shock.dates, col = "blue",  lty = 1, lwd = 2.5)
  if (legend) legend("topleft", legend = colnames(x), col = 1:ncol(m), lty = 1, cex = 0.6)
}
plot_series(prices_clean, "Cleaned Adjusted Prices", "Price", legend = TRUE)
plot_series(log_returns,  "log-returns",             "Log-return")
par(mfrow = c(1, 1))

X.daily <- as.matrix(returns_in)
MardiaTest(X.daily)

JBTest <- function(dataMatrix)
{
  pvals <- rep(NA, ncol(dataMatrix))
  for (i in 1:ncol(dataMatrix))
  {
    pvals[i] <- jarque.bera.test(dataMatrix[, i])$p.value
  }
  names(pvals) <- colnames(dataMatrix)
  return(pvals)
}
JBTest(X.daily)

jb.stats <- apply(X.daily, 2, function(x) jarque.bera.test(x)$statistic)
worst <- names(which.max(jb.stats))
qqnorm(X.daily[, worst], main = paste("Normal QQ-plot:", worst, "daily log-returns"))
qqline(X.daily[, worst], col = "red", lwd = 2)

moments <- t(apply(X.daily, 2, function(x) {
  s <- mean((x - mean(x))^3) / (mean((x - mean(x))^2))^(3/2)
  k <- mean((x - mean(x))^4) / (mean((x - mean(x))^2))^2
  c(skewness = s, kurtosis = k)
}))
round(moments, 2)

X.weekly    <- as.matrix(apply.weekly(X.daily, colSums))
X.monthly   <- as.matrix(apply.monthly(X.daily, colSums))
X.quarterly <- as.matrix(apply.quarterly(X.daily, colSums))

round(MardiaTest(X.weekly), 4)
round(MardiaTest(X.monthly), 4)
round(MardiaTest(X.quarterly), 4)

jb.results <- rbind(
  daily     = JBTest(X.daily),
  weekly    = JBTest(X.weekly),
  monthly   = JBTest(X.monthly),
  quarterly = JBTest(X.quarterly))
round(jb.results, 4)
prop.reject <- apply(jb.results, 1, function(p) mean(p < 0.05))
prop.reject
#--------  hier boven compiles goed -------:)----
consumer   <- c("MC.PA", "OR.PA", "ADS.DE")
technology <- c("SAP.DE", "ASML.AS")
energy     <- c("TTE.PA", "ENI.MI")
banking    <- c("BNP.PA", "ISP.MI")
across     <- c("MC.PA", "SAP.DE", "SIE.DE", "TTE.PA", "BNP.PA")

panel.cor <- function(x, y, ...) {
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.5, sprintf("%.2f", cor(x, y, method = "spearman")), cex = 1.5)
}

pairs(as.data.frame(X.daily[, consumer]),
      upper.panel = points, lower.panel = panel.cor,
      main = "Within sector: Consumer", pch = 20, cex = 0.4)
pairs(as.data.frame(X.daily[, technology]),
      upper.panel = points, lower.panel = panel.cor,
      main = "Within sector: Technology", pch = 20, cex = 0.4)
pairs(as.data.frame(X.daily[, energy]),
      upper.panel = points, lower.panel = panel.cor,
      main = "Within sector: Energy", pch = 20, cex = 0.4)
pairs(as.data.frame(X.daily[, banking]),
      upper.panel = points, lower.panel = panel.cor,
      main = "Within sector: Banking", pch = 20, cex = 0.4)
pairs(as.data.frame(X.daily[, across]),
      upper.panel = points, lower.panel = panel.cor,
      main = "Across sectors: one stock per sector", pch = 20, cex = 0.4)

cor.mat   <- cor(X.daily, method = "spearman")
sector.of <- setNames(sector_map$sector, sector_map$ticker)
all.pairs   <- t(combn(colnames(X.daily), 2))
pair.cors   <- cor.mat[all.pairs]
same.sector <- sector.of[all.pairs[, 1]] == sector.of[all.pairs[, 2]]
cat("Average within-sector correlation:", round(mean(pair.cors[same.sector]), 2), "\n")
cat("Average across-sector correlation:", round(mean(pair.cors[!same.sector]), 2), "\n")
#____
initial_prices <- as.numeric(prices_clean[1, ])
lambda <- 1000 / initial_prices
names(lambda) <- colnames(prices_clean)
V <- xts(coredata(prices_clean) %*% lambda, order.by = index(prices_clean))
colnames(V) <- "Value"
L <- -diff(V)[-1, ]
colnames(L) <- "loss"

d <- index(V); v <- as.numeric(V)
plot(d, v, type = "l", xlab = "", ylab = "Value (EUR)",
     main = "Portfolio value, 1000 EUR per stock, never rebalanced", xaxt = "n")
axis.Date(1, at = seq(min(d), max(d), by = "year"), format = "%Y")
abline(v = split.date,  col = "black", lty = 2, lwd = 2)
abline(v = shock.dates, col = "blue",  lty = 1, lwd = 2)
abline(h = 10000, col = "grey60", lty = 3)

L_in  <- L[paste(start_date, end_date_in_sample, sep = "/")]
L_out <- L[paste(start_date_out_sample, end_date, sep = "/")]
L.in  <- as.numeric(L_in)
L.out <- as.numeric(L_out)
plot(L, main = "Daily portfolio loss (EUR)")

VaR_emp <- function(a) quantile(L.in, a, names = FALSE)
ES_emp  <- function(a) mean(L.in[L.in >= VaR_emp(a)])

mu.n  <- mean(L.in)
sig.n <- sd(L.in)
VaR_norm <- function(a) mu.n + sig.n * qnorm(a)
ES_norm  <- function(a) mu.n + sig.n * dnorm(qnorm(a)) / (1 - a)
#_____
tfit <- fit.st(L.in)
nu.t  <- tfit$par.ests["nu"]
mu.t  <- tfit$par.ests["mu"]
sig.t <- tfit$par.ests["sigma"]
VaR_t <- function(a) mu.t + sig.t * qt(a, df = nu.t)
ES_t  <- function(a) {
  ta <- qt(a, df = nu.t)
  ES.std <- dt(ta, df = nu.t) / (1 - a) * (nu.t + ta^2) / (nu.t - 1)
  mu.t + sig.t * ES.std
}

alphas <- c(0.95, 0.99)
risk.table <- data.frame(
  alpha    = alphas,
  VaR_emp  = sapply(alphas, VaR_emp),
  VaR_norm = sapply(alphas, VaR_norm),
  VaR_t    = sapply(alphas, VaR_t),
  ES_emp   = sapply(alphas, ES_emp),
  ES_norm  = sapply(alphas, ES_norm),
  ES_t     = sapply(alphas, ES_t))
round(risk.table, 2)

T.out <- length(L.out)
count_viol <- function(VaR.value) sum(L.out >= VaR.value)

backtest <- data.frame(
  alpha    = alphas,
  expected = ceiling((1 - alphas) * T.out),
  emp      = sapply(sapply(alphas, VaR_emp),  count_viol),
  norm     = sapply(sapply(alphas, VaR_norm), count_viol),
  t        = sapply(sapply(alphas, VaR_t),    count_viol))
backtest

kupiec <- function(I, T, a) {
  phat <- I / T
  if (I == 0) {
    LR <- 2 * ((T - I) * log((1 - phat) / a))
  } else {
    LR <- 2 * (I * log(phat / (1 - a)) + (T - I) * log((1 - phat) / a))
  }
  pval <- 1 - pchisq(LR, df = 1)
  c(violations = I, expected = (1 - a) * T, LR = LR, p.value = pval)
}

models <- list(emp = VaR_emp, norm = VaR_norm, t = VaR_t)

for (a in alphas) {
  cat("\n---- alpha =", a, "----\n")
  for (m in names(models)) {
    I <- count_viol(models[[m]](a))
    res <- kupiec(I, T.out, a)
    cat(sprintf("%-5s violations = %2d (expected %.1f)  LR = %6.2f  p = %.4f\n",
                m, res["violations"], res["expected"], res["LR"], res["p.value"]))
  }
}

for (a in alphas) {
  cat("\n---- alpha =", a, "----\n")
  for (m in names(models)) {
    viol   <- as.numeric(L.out >= models[[m]](a))
    tab    <- table(prev = viol[-length(viol)], now = viol[-1])
    p.test <- chisq.test(tab, correct = FALSE)
    cat(sprintf("%-5s P(viol | viol yesterday) = %.3f  P(viol | calm yesterday) = %.3f  p = %.4f\n",
                m, tab["1","1"] / sum(tab["1", ]), tab["0","1"] / sum(tab["0", ]),
                p.test$p.value))
  }
}

plot(index(L_out), L.out, type = "h", col = "grey40",
     xlab = "", ylab = "Daily loss (EUR)",
     main = "Backtest 2022-2024: losses vs in-sample 99% VaR")
abline(h = VaR_emp(0.99),  col = "green", lwd = 2)
abline(h = VaR_norm(0.99), col = "red",   lwd = 2, lty = 2)
abline(h = VaR_t(0.99),    col = "blue",  lwd = 2, lty = 3)
legend("topright", c("empirical", "normal", "Student-t"),
       col = c("green", "red", "blue"), lty = 1:3, lwd = 2, cex = 0.8)

X <- as.matrix(returns_in)
X <- X[apply(X != 0, 1, all), ]
n <- nrow(X)
n
copulaX <- apply(X, 2, edf, adjust = 1)

pairs(as.data.frame(copulaX), upper.panel = points, lower.panel = panel.cor,
      xlim = c(0, 1), ylim = c(0, 1),
      main = "Pseudo copula data: all 10 stocks pairwise",
      pch = 20, cex = 0.2)

cor_s <- cor(copulaX, method = "spearman")
pairs_all <- t(combn(colnames(copulaX), 2))
ranking <- data.frame(pair = paste(pairs_all[,1], pairs_all[,2], sep = "-"),
                      spearman = round(cor_s[pairs_all], 2))
head(ranking[order(-ranking$spearman), ], 5)

pair <- c("BNP.PA", "ISP.MI")
U <- copulaX[, pair]

fitGauss <- fit.gausscopula(U)
fitT     <- fit.tcopula(U)
fitGumb  <- fit.AC(U, "gumbel")
fitClay  <- fit.AC(U, "clayton")

fitGauss$P[1,2]
fitT$P[1,2]
fitT$nu
fitGumb$theta
fitClay$theta

ll <- c(fitGauss$ll.max, fitT$ll.max, fitGumb$ll.max, fitClay$ll.max)
k  <- c(1, 2, 1, 1)
comparison <- data.frame(
  copula = c("Gauss", "StudT", "Gumbel", "Clayton"),
  logLik = round(ll, 1),
  AIC    = round(-2 * ll + 2 * k, 1),
  BIC    = round(-2 * ll + k * log(n), 1))
comparison

pairs_list <- list(c("MC.PA","OR.PA"),
                   c("SAP.DE","ASML.AS"),
                   c("TTE.PA","ENI.MI"),
                   c("BNP.PA","ISP.MI"),
                   c("MC.PA","SAP.DE"))
best <- data.frame(pair = character(), best_copula = character(), t_nu = numeric())
for (p in pairs_list) {
  Up <- copulaX[, p]
  fG <- fit.gausscopula(Up); fT <- fit.tcopula(Up)
  fGu <- fit.AC(Up, "gumbel"); fCl <- fit.AC(Up, "clayton")
  aic <- -2 * c(fG$ll.max, fT$ll.max, fGu$ll.max, fCl$ll.max) + 2 * c(1,2,1,1)
  best <- rbind(best, data.frame(
    pair = paste(p, collapse = "-"),
    best_copula = c("Gauss","StudT","Gumbel","Clayton")[which.min(aic)],
    t_nu = round(fT$nu, 1)))
}
best

fitGauss10 <- fit.gausscopula(copulaX)
fitT10     <- fit.tcopula(copulaX)
fitT10$nu

d <- ncol(copulaX)
k10 <- c(d*(d-1)/2, d*(d-1)/2 + 1)
ll10 <- c(fitGauss10$ll.max, fitT10$ll.max)
data.frame(copula = c("Gauss", "StudT"),
           logLik = round(ll10, 1),
           AIC = round(-2*ll10 + 2*k10, 1),
           BIC = round(-2*ll10 + k10*log(n), 1))

rho <- fitT$P[1,2]; nu <- fitT$nu
lambda_t <- 2 * pt(-sqrt((nu + 1) * (1 - rho) / (1 + rho)), df = nu + 1)
lambda_t

for (q in c(0.10, 0.05, 0.025)) {
  emp <- sum(U[,1] < q & U[,2] < q) / sum(U[,1] < q)
  cat("q =", q, ": empirical P(joint crash | one crashes) =", round(emp, 3), "\n")
}