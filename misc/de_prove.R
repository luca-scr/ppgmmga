library(ppgmmga)
library(mclust)
library(ggplot2)
theme_set(theme_bw())
# theme_update(plot.title = element_text(hjust = 0.5))

# gaControl("useRcpp" = FALSE)

# library(RcppDE)
# ?DEoptim

# Waveform data -----------------------------------------------------------

library(mlbench)
set.seed(20180124)
x <- mlbench.waveform(400)
X <- x$x
Class <- factor(x$classes)

gmm <- densityMclust(scale(X, center = TRUE, scale = FALSE))

PPGMMGA1 <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class, drawAxis = FALSE) + ggtitle("UT[GA]")

PPGMMDE <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1, gatype = "de", opt = list(optim = FALSE))

PPGMMDE <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1, gatype = "de", opt = list(optim = FALSE, maxiter = 100, keepBest = TRUE))
GA::plot(PPGMMDE$GA)
summary(PPGMMDE$GA)
summary(PPGMMDE, check = TRUE)
plot(PPGMMDE, Class, drawAxis = FALSE) + ggtitle("UT[DE]")

which(diff(PPGMMDE$GA@summary[,1]) < 0)
PPGMMDE$GA@bestSol

# Crabs data -----------------------------------------------------------

data(crabs, package = "MASS")
X <- crabs[, 4:8]
Class <- as.factor(with(crabs, paste(sp, sex, sep = "|")))

X <- scale(X, center = TRUE, scale = FALSE)
gmm <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, d = 2, gmm = gmm, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class) + ggtitle("UT")

PPGMMDE <- ppgmmga(data = X, d = 2, gmm = gmm, approx = "UT", seed = 1, gatype = "de")
summary(PPGMMDE, check = TRUE)
plot(PPGMMDE, Class) + ggtitle("UT[DE]")


# Coffee data ----------------------------------------------------------

data("coffee", package = "pgmm")
X <- coffee[,-(1:2)]
names(X)[8] <- c("Caffeine")
Class <- factor(coffee$Variety, levels = 1:2, labels = c("Arabica", "Robusta"))

X <- scale(X, center = TRUE)
gmm <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, d = 2, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class) + ggtitle("UT")

PPGMMDE <- ppgmmga(data = X, d = 2, approx = "UT", seed = 1, gatype = "de")
PPGMMDE <- ppgmmga(data = X, d = 2, approx = "UT", seed = 1, gatype = "de",
                   opt = list(optim = TRUE, run = 200, pmutation = 0.1, elitism = 1))
plot(PPGMMDE$GA)
summary(PPGMMDE, check = TRUE)
plot(PPGMMDE, Class) + ggtitle("UT[DE]")


# AIS data -------------------------------------------------------------

data(ais, package = "dr")
X <- ais[,2:12]
Class <- factor(ifelse(ais$Sex == 0, "M", "F"))

X <- scale(X, center = TRUE)
gmm <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class) + ggtitle("UT")

PPGMMDE <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1, gatype = "de")
PPGMMDE <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1, gatype = "de",
                   opt = list(optim = TRUE, run = 200, pmutation = 0.1, elitism = 1))
plot(PPGMMDE$GA)
summary(PPGMMDE, check = TRUE)
plot(PPGMMDE, Class) + ggtitle("UT[DE]")


# Banknote authentication data -----------------------------------------

banknote <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/00267/data_banknote_authentication.txt", header = FALSE)
colnames(banknote) <- c("WavVar", "WavSkw", "WavKur", "Entropy", "Class")
X <- banknote[,1:4]
Class <- as.factor(ifelse(banknote$Class == 0, "genuine", "forged"))

X <- scale(X, center = TRUE)
gmm <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class) + ggtitle("UT")

PPGMMDE <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1, gatype = "de")
PPGMMDE <- ppgmmga(data = X, gmm = gmm, d = 2, approx = "UT", seed = 1, gatype = "de",
                   opt = list(optim = FALSE, run = 200, pmutation = 0.1, elitism = 1))
plot(PPGMMDE$GA)
summary(PPGMMDE, check = TRUE)
plot(PPGMMDE, Class) + ggtitle("UT[DE]")



