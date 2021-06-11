# Install packages.
# install.packages("tidyverse")

# Load libraries.
library("ggplot2")

# Set directory.
# setwd("C:/Users/Stefan van Berkum/git/HAPSA/HAPSA/Results/")
setwd("C:/Users/Stefan van Berkum/Google Drive/Studie/Thesis/Back-up/Results Old/")
plot_out <- "C:/Users/Stefan van Berkum/Google Drive/Studie/Thesis/Back-up/Graphs Old/"

# Read files.
apsa <- read.table("APSA/CSV/30_240_scores.csv", header = TRUE, sep = ",")

gamm <- c(
  0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0,
  115.0, 120.0, 125.0
)

results <- apsa

for (i in 2:length(gamm)) {
  varName <- paste("vis", gamm[i], sep = "_")
  fileName <- paste("VIS/CSV/30_240_", gamm[i], ".0_scores.csv", sep = "")
  table <- assign(varName, read.table(fileName, header = TRUE, sep = ","))
  results <- rbind(results, table)
}
results <- cbind(gamm, results)
names(results) <- cbind("Gamma", "Profit", "SHS", "SVHS")

# Calculate deltas.
fdVar <- function(x, var) {
  (x - apsa[var]) / apsa[var] * 100
}
dProfit <- data.frame(unlist(sapply(results[, 2], fdVar, var = "Profit")))
dSHS <- data.frame(unlist(sapply(results[, 3], fdVar, var = "SHS")))
dSVHS <- data.frame(unlist(sapply(results[, 4], fdVar, var = "SVHS")))
dProfit <- cbind(gamm, dProfit, "Profit")
names(dProfit) <- cbind("Gamma", "Change", "Score")
dSHS <- cbind(gamm, dSHS, "SHS")
names(dSHS) <- cbind("Gamma", "Change", "Score")
dSVHS <- cbind(gamm, dSVHS, "SVHS")
names(dSVHS) <- cbind("Gamma", "Change", "Score")
dResults <- rbind(dProfit, dSHS, dSVHS)

# Calculate ratios dSHS/dProfit and dSVHS/dProfit.
dProfitVec <- data.frame(unlist(sapply(results[, 2], fdVar, var = "Profit")))
dProfitVec <- data.frame(dProfitVec[-1,])
dSHSVec <- data.frame(unlist(sapply(results[, 3], fdVar, var = "SHS")))
dSHSVec <- data.frame(dSHSVec[-1,])
dSVHSVec <- data.frame(unlist(sapply(results[, 4], fdVar, var = "SVHS")))
dSVHSVec <- data.frame(dSVHSVec[-1,])
dSHS_dProfit <- dSHSVec / abs(dProfitVec)
SHS_ratios <- cbind(data.frame(gamm)[-1,], dSHS_dProfit, "dSHS / |dProfit|")
names(SHS_ratios) <- cbind("Gamma", "Ratio", "Score")
dSVHS_dProfit <- dSVHSVec / abs(dProfitVec)
SVHS_ratios <- cbind(data.frame(gamm)[-1,], dSVHS_dProfit, "dSVHS / |dProfit|")
names(SVHS_ratios) <- cbind("Gamma", "Ratio", "Score")
ratios <- rbind(SHS_ratios, SVHS_ratios)

# Calculate first differences delta(dProfit), delta(dSHS), and delta(dSVHS)
# using an x-observation moving window.
window_size <- 5
first <- window_size + 1
partialGamma <- data.frame(gamm)[first:length(gamm),]
mw <- mean(dProfit[1:window_size, 2])
fd_dProfit <- data.frame(dProfit[first, 2] - mw)
names(fd_dProfit) <- "Change"
for (i in (first + 1):length(dProfit[,2])) {
  window_start <- i - window_size
  window_end <- i - 1
  mw <- mean(dProfit[window_start:window_end, 2])
  fd <- data.frame(dProfit[i, 2] - mw)
  names(fd) <- "Change"
  fd_dProfit <- rbind(fd_dProfit, fd)
}
fd_dProfit <- cbind(partialGamma, abs(fd_dProfit), "Profit")
names(fd_dProfit) <- cbind("Gamma", "Change", "Score")

mw <- mean(dSHS[1:window_size, 2])
fd_dSHS <- data.frame(dSHS[first, 2] - mw)
names(fd_dSHS) <- "Change"
for (i in (first + 1):length(dSHS[,2])) {
  window_start <- i - window_size
  window_end <- i - 1
  mw <- mean(dSHS[window_start:window_end, 2])
  fd <- data.frame(dSHS[i, 2] - mw)
  names(fd) <- "Change"
  fd_dSHS <- rbind(fd_dSHS, fd)
}
fd_dSHS <- cbind(partialGamma, abs(fd_dSHS), "SHS")
names(fd_dSHS) <- cbind("Gamma", "Change", "Score")

mw <- mean(dSVHS[1:window_size, 2])
fd_dSVHS <- data.frame(dSVHS[first, 2] - mw)
names(fd_dSVHS) <- "Change"
for (i in (first + 1):length(dSVHS[,2])) {
  window_start <- i - window_size
  window_end <- i - 1
  mw <- mean(dSVHS[window_start:window_end, 2])
  fd <- data.frame(dSVHS[i, 2] - mw)
  names(fd) <- "Change"
  fd_dSVHS <- rbind(fd_dSVHS, fd)
}
fd_dSVHS <- cbind(partialGamma, fd_dSVHS, "SVHS")
names(fd_dSVHS) <- cbind("Gamma", "Change", "Score")

fd_dResults = rbind(fd_dProfit, fd_dSHS, fd_dSVHS)

# Plot results.
# Profit.
ggplot(data = results, aes(x = Gamma, y = Profit)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = unlist(apsa["Profit"]))) +
  xlab(expression(gamma))
out_path = paste(plot_out, "Profit.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# SHS.
ggplot(data = results, aes(x = Gamma, y = SHS)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = unlist(apsa["SHS"]))) +
  xlab(expression(gamma))
out_path = paste(plot_out, "SHS.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# SVHS.
ggplot(data = results, aes(x = Gamma, y = SVHS)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = unlist(apsa["SVHS"]))) +
  xlab(expression(gamma))
out_path = paste(plot_out, "SVHS.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# dResults.
ggplot(data = dResults, aes(x = Gamma, y = Change, color = Score, shape = Score)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = 0)) +
  xlab(expression(gamma)) +
  ylab("Change (%)")
out_path = paste(plot_out, "dResults.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# dSHS/dProfit.
ggplot(data = SHS_ratios, aes(x = Gamma, y = Ratio)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = 1)) +
  xlab(expression(gamma))
out_path = paste(plot_out, "dSHS_dProfit.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# dSVHS/dProfit.
ggplot(data = SVHS_ratios, aes(x = Gamma, y = Ratio)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = 1)) +
  xlab(expression(gamma))
out_path = paste(plot_out, "dSVHS_dProfit.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# dSHS/dProfit and dSVHS/dProfit.
ggplot(data = ratios, aes(x = Gamma, y = Ratio, color = Score, shape = Score)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  geom_hline(aes(yintercept = 1)) +
  xlab(expression(gamma))
out_path = paste(plot_out, "dResults_dProfit.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)

# First differences delta(dProfit), delta(dSHS), and delta(dProfit).
ggplot(data = fd_dResults, aes(x = Gamma, y = Change, color = Score, shape = Score)) +
  geom_point() +
  geom_smooth(method = loess, formula = y ~ x) +
  xlab(expression(gamma)) +
  ylab(expression(paste(Delta, "|Change|")))
out_path = paste(plot_out, "delta_dResults.png")
ggsave(out_path, width = 150, height = 100, units = "mm", dpi = 600)
