#' ---
#' title: 'R code accompanying the paper "Projection pursuit based on Gaussian mixtures and evolutionary algorithms", JCGS, 2019'
#' author: "Luca Scrucca and Alessio Serafini"
#' date: "15 Oct 2018"
#' output:
#'   pdf_document:
#'     highlight: null
#'     fig_width: 6
#'     fig_height: 5
#' ---

# If needed please install the following packages from CRAN:
#
# install.packages(c("ppgmmga", "gridExtra", "rmarkdown", "mlbench", "gtable",
#                    "pgmm", "dr", "fastICA"), dependencies = TRUE)
#
# and the following packages from Bioconductor:
# 
# if(!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("Biobase")
# BiocManager::install("multtest")

# To reproduce the results of the paper please use:
# - R ver. 3.5.3 (2019-03-11) 
# - ppgmmga ver. 1.2 (2019-07-08)
# In R ver. >= 3.6 the following command is needed for backward compatibility:
RNGkind(sample.kind = "Rounding")

# To compile the full report with results use:
# rmarkdown::render("code.R")

library(ppgmmga)
library(mclust)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))
source(system.file("scripts/jcgs", "mc_negent_other_methods.R", package = "ppgmmga"))
source(system.file("scripts/jcgs", "misc.R", package = "ppgmmga"))


# Waveform data -----------------------------------------------------------

library(mlbench)
set.seed(20180124)
x <- mlbench.waveform(400)
X <- x$x
Class <- factor(x$classes)

X <- scale(X, center = TRUE, scale = FALSE)
GMM <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, d = 2, GMM = GMM, scale = FALSE, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class, drawAxis = FALSE) + ggtitle("PPGMMGA[UT]")  

PPGMMGA2 <- ppgmmga(data = X, d = 2, GMM = GMM, scale = FALSE, approx = "VAR", seed = 2)
summary(PPGMMGA2, check = TRUE)
plot(PPGMMGA2, Class, drawAxis = FALSE) + ggtitle("PPGMMGA[VAR]") 

PPGMMGA3 <- ppgmmga(data = X, d = 2, GMM = GMM, scale = FALSE, approx = "SOTE", seed = 3)
summary(PPGMMGA3, check = TRUE)
plot(PPGMMGA3, Class, drawAxis = FALSE) + ggtitle("PPGMMGA[SOTE]") 

PCA <- NegentropyPCA(PPGMMGA1)
PCA[c("Negentropy", "se")]
PPGMMPCA <- PPGMMGA1; PPGMMPCA$approx <- "PCA"
PPGMMPCA$basis <- PCA$basis
PPGMMPCA$Z <- PCA$Z
plot(PPGMMPCA, Class, drawAxis = FALSE) + ggtitle("PCA")

ICA <- NegentropyFASTICA(PPGMMGA1)
ICA[c("Negentropy", "se")]
# trick for plotting
PPGMMICA <- PPGMMGA1; PPGMMICA$approx <- "ICA"
PPGMMICA$basis <- ICA$basis
PPGMMICA$Z <- ICA$Z
plot(PPGMMICA, Class, drawAxis = FALSE) + ggtitle("ICA")

ppgmmga.table.angle(PPGMMGA1, PPGMMGA2, PPGMMGA3, PPGMMPCA, PPGMMICA)

# Crabs data -----------------------------------------------------------

data(crabs, package = "MASS")
X <- crabs[, 4:8]
Class <- as.factor(with(crabs, paste(sp, sex, sep = "|")))

X <- scale(X, center = TRUE, scale = TRUE)
GMM <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, d = 2, GMM = GMM, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class) + ggtitle("PPGMMGA[UT]")  

PPGMMGA2 <- ppgmmga(data = X, d = 2, GMM = GMM, approx = "VAR", seed = 2)
summary(PPGMMGA2, check = TRUE)
plot(PPGMMGA2, Class) + ggtitle("PPGMMGA[VAR]") 

PPGMMGA3 <- ppgmmga(data = X, d = 2, GMM = GMM, approx = "SOTE", seed = 3)
summary(PPGMMGA3, check = TRUE)
plot(PPGMMGA3, Class) + ggtitle("PPGMMGA[SOTE]") 

PCA <- NegentropyPCA(PPGMMGA1)
PCA[c("Negentropy", "se")]
PPGMMPCA <- PPGMMGA1; PPGMMPCA$approx <- "PCA"
PPGMMPCA$basis <- PCA$basis
PPGMMPCA$Z <- PCA$Z
plot(PPGMMPCA, Class) + ggtitle("PCA")

ICA <- NegentropyFASTICA(PPGMMGA1)
ICA[c("Negentropy", "se")]
# trick for plotting
PPGMMICA <- PPGMMGA1; PPGMMICA$approx <- "ICA"
PPGMMICA$basis <- ICA$basis
PPGMMICA$Z <- ICA$Z
plot(PPGMMICA, Class, drawAxis = FALSE) + ggtitle("ICA")

ppgmmga.table.angle(PPGMMGA1, PPGMMGA2, PPGMMGA3, PPGMMPCA, PPGMMICA)


# Coffee data ----------------------------------------------------------

data("coffee", package = "pgmm")
X <- coffee[,-(1:2)]
names(X)[8] <- c("Caffeine")
Class <- factor(coffee$Variety, levels = 1:2, labels = c("Arabica", "Robusta"))

X <- scale(X, center = TRUE, scale = TRUE)
GMM <- densityMclust(X)

PPGMMGA1 <- ppgmmga(data = X, d = 1, GMM = GMM, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class, bins = 9) + ggtitle("PPGMMGA[UT]")  

PPGMMGA2 <- ppgmmga(data = X, d = 1, GMM = GMM, approx = "VAR", seed = 2)
summary(PPGMMGA2, check = TRUE)
plot(PPGMMGA2, Class, bins = 9) + ggtitle("PPGMMGA[VAR]") 

PPGMMGA3 <- ppgmmga(data = X, d = 1, GMM = GMM, approx = "SOTE", seed = 3)
summary(PPGMMGA3, check = TRUE)
plot(PPGMMGA3, Class, bins = 9) + ggtitle("PPGMMGA[SOTE]") 

PCA <- NegentropyPCA(PPGMMGA1)
PCA[c("Negentropy", "se")]
PPGMMPCA <- PPGMMGA1; PPGMMPCA$approx <- "PCA"
PPGMMPCA$basis <- PCA$basis
PPGMMPCA$Z <- PCA$Z
plot(PPGMMPCA, Class, nbins = 9) + ggtitle("PCA")

ICA <- NegentropyFASTICA(PPGMMGA1)
ICA[c("Negentropy", "se")]
# trick for plotting
PPGMMICA <- PPGMMGA1; PPGMMICA$approx <- "ICA"
PPGMMICA$basis <- ICA$basis
PPGMMICA$Z <- ICA$Z
plot(PPGMMICA, Class, bins = 9) + ggtitle("ICA")

ppgmmga.table.angle(PPGMMGA1, PPGMMGA2, PPGMMGA3, PPGMMPCA, PPGMMICA)

df <- data.frame(variable = rownames(PPGMMGA1$basis), coefs = PPGMMGA1$basis[,1])
plot1 <- ggplot(df, aes(x = variable)) + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(ymin = ifelse(coefs < 0, coefs, 0),
                     ymax = ifelse(coefs > 0, coefs, 0)),
                 lwd = 1, position = position_dodge(width = 1/2)) + 
  xlab("") + ylab("PP1 coefficients") +
  coord_flip() +
  theme_bw()
df <- cbind(data.frame(X, check.names = FALSE)[c("Caffeine", "Fat")], Class)
plot2 <- ggplot(df, aes(Class, Caffeine, fill = Class, color = Class)) +
  geom_boxplot(outlier.shape = 19, alpha = 1/2) + 
  scale_fill_tableau("Classic 10") +
  scale_colour_tableau("Classic 10")
plot3 <- ggplot(df, aes(Class, Fat, fill = Class, color = Class)) +
  geom_boxplot(outlier.shape = 19, alpha = 1/2) + 
  scale_fill_tableau("Classic 10") +
  scale_colour_tableau("Classic 10")
plots <- grid.arrange(plot1, 
                      multiple_ggplot_sharedLegend(plot2, plot3, nrow = 1, 
                                                   position = "bottom"),
                      nrow = 1, widths = c(1,1))


# AIS data -------------------------------------------------------------

data(ais, package = "dr")
X <- ais[,2:12]
Class <- factor(ifelse(ais$Sex == 0, "M", "F"))

X <- scale(X, center = TRUE, scale = TRUE)
GMM <- densityMclust(X)

# d = 1
PPGMMGA1 <- ppgmmga(data = X, d = 1, GMM = GMM, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)

PPGMMGA2 <- ppgmmga(data = X, d = 1, GMM = GMM, approx = "VAR", seed = 2)
summary(PPGMMGA2, check = TRUE)

PPGMMGA3 <- ppgmmga(data = X, d = 1, GMM = GMM, approx = "SOTE", seed = 3)
summary(PPGMMGA3, check = TRUE)

PCA <- NegentropyPCA(PPGMMGA1)
PCA[c("Negentropy", "se")]
PPGMMPCA <- PPGMMGA1; PPGMMPCA$approx <- "PCA"
PPGMMPCA$basis <- PCA$basis
PPGMMPCA$Z <- PCA$Z

ppgmmga.table.angle(PPGMMGA1, PPGMMGA2, PPGMMGA3, PPGMMPCA)

# d = 2
PPGMMGA1 <- ppgmmga(data = X, d = 2, GMM = GMM, approx = "UT", seed = 1)
summary(PPGMMGA1, check = TRUE)
plot(PPGMMGA1, Class) + ggtitle("PPGMMGA[UT]") 

PPGMMGA2 <- ppgmmga(data = X, d = 2, GMM = GMM, approx = "VAR", seed = 2)
summary(PPGMMGA2, check = TRUE)
plot(PPGMMGA2, Class) + ggtitle("PPGMMGA[VAR]") 

PPGMMGA3 <- ppgmmga(data = X, d = 2, GMM = GMM, approx = "SOTE", seed = 3)
summary(PPGMMGA3, check = TRUE)
plot(PPGMMGA3, Class) + ggtitle("PPGMMGA[SOTE]") 

PCA <- NegentropyPCA(PPGMMGA1)
PCA[c("Negentropy", "se")]
PPGMMPCA <- PPGMMGA1; PPGMMPCA$approx <- "PCA"
PPGMMPCA$basis <- PCA$basis
PPGMMPCA$Z <- PCA$Z
plot(PPGMMPCA, Class) + ggtitle("PCA")

ppgmmga.table.angle(PPGMMGA1, PPGMMGA2, PPGMMGA3, PPGMMPCA)


# Leukemia data ---------------------------------------------------------

data(golub, package = "multtest")
X <- t(as.matrix(golub))
X <- scale(X, center = TRUE, scale = FALSE)
dim(X)
Class <- factor(golub.cl, levels = 0:1, labels = c("ALL", "AML"))
table(Class)

GMM <- densityMclust(X, modelNames = c("EII", "VII", "EEI", "EVI", "VEI", "VVI"))
summary(GMM)
m <- GMM$parameters$mean
s <- sqrt(apply(GMM$parameters$variance$sigma, 3, diag))
Signal <- (m[,1]-m[,2])
Signal2Noise <- (m[,1]-m[,2])/(s[,1]+s[,2])
expr <- (abs(Signal) > 1 & abs(Signal2Noise) > 1)
plot(Signal, abs(Signal2Noise), 
     cex = 0.5, pch = ifelse(expr, 19, 1))
abline(h = c(-1,1), v = c(-1,1), lty = 2)
g <- which(expr)

GMM2 <- densityMclust(X[,g], modelNames = c("EII", "VII", "EEI", "EVI", "VEI", "VVI"))
PPGMMGA <- ppgmmga(data = X[,g], d = 2, approx = "UT", seed = 1,
                   GMM = GMM2, scale = FALSE, 
                   options = ppgmmga.options(maxiter = 2000))
summary(PPGMMGA, check = TRUE)
plot(PPGMMGA, Class)

df <- data.frame(Signal, Signal2Noise, expr)
plot1 <- ggplot(df, aes(x = Signal, y = abs(Signal2Noise), color = expr)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  scale_x_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3,3)) +
  xlab("Signal (difference of means)") + 
  ylim(c(0, 2)) + 
  ylab("Abs signal to noise ratio") +
  theme(legend.position = "none")

plot2 <- 
plot(PPGMMGA, Class, drawAxis = FALSE) + 
  scale_x_continuous(expand = expand_scale(c(0.5, 0.1))) +
  scale_y_continuous(expand = expand_scale(c(0.2, 0.1))) +
  guides(col = guide_legend(title = "Leukemia type"),
         pch = guide_legend(title = "Leukemia type")) +
  theme(legend.position = c(0.02, 0.02), 
        legend.direction = "horizontal",
        legend.justification = c(0,0),
        legend.margin = margin(2, 2, 2, 2),
        legend.box.background = element_rect(colour = "black"))

plots <- grid.arrange(plot1, plot2, nrow = 1, ncol = 2)

# ----

# Restore default sample's behavior:
RNGkind(sample.kind = "default")
