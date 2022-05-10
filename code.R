# jonashaslbeck@gmail.com; May 2, 2022

# ----------------------------------------------------------
# -------- What is this? -----------------------------------
# ----------------------------------------------------------

# Do the analysis for my session for NSMD workshop in Rowland data


# ----------------------------------------------------------
# -------- Load Packages -----------------------------------
# ----------------------------------------------------------

library(plyr)
library(RColorBrewer)
source("fitVAR.R")
library(moments)
library(mgm)
library(qgraph)
library(mlVAR)


# ----------------------------------------------------------
# -------- Aux functions -----------------------------------
# ----------------------------------------------------------


# Aux function for plotting:
PL <- function(tex, cex = 1.5, srt=0, x=0.45, y=0.5) {
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(x, y, labels = tex, srt=srt, cex=cex, adj=0.25)
}


# ----------------------------------------------------------
# -------- Load Data ---------------------------------------
# ----------------------------------------------------------

data <- readRDS("ESMdata.RDS")
head(data)

# ----------------------------------------------------------
# -------- Missing data ------------------------------------
# ----------------------------------------------------------

# Average missingness across people
out <- ddply(data, .(subj_id), function(x) {
  mean(!is.na(x$happy))
})
hist(out$V1, xlim=c(0,1), xlab="Propotion Missing",
     breaks=seq(0, 1, length=30), main="")


# Average missingness across people
cols <- brewer.pal(5, "Set1")
plot.new()
plot.window(xlim=c(1,240), ylim=c(1,240))
axis(1, c(1, 50 ,100, 150, 200, 240))
axis(2, las=2, c(1, 50 ,100, 150, 200, 240))
title(xlab="Time steps", ylab="Number of measurements")
title(main="Cumulative measurements for 5 subjects", font.main = 1)
lines(1:240, col="grey", lwd=2)
u_subj <- unique(data$subj_id)
for(i in 1:5) {
  data_i <- data[data$subj_id==u_subj[i], ]
  miss_seq <- !is.na(data_i$excited)
  miss_seq_cs <- cumsum(miss_seq)
  lines(miss_seq_cs, col=cols[i], lwd=2)
}
legend("topleft", legend=c("No Missing", paste0("Subject ", 1:5)),
       col=c("grey", cols), lwd=rep(2, 6), bty="n")

# Length of Missing Sequences
l_ints <- list()
for(i in 1:125) {
  data_i <- data[data$subj_id==u_subj[i], ]
  v_NAseq <- c()
  c1 <- 0

  for(j in 1:239) {

    if(is.na(data_i$happy[j])) {
      c1 <- c1 + 1
    } else {
      v_NAseq <- c(v_NAseq, c1)
      c1 <- 0
    }

  } # end for: j

  l_ints[[i]] <- v_NAseq

} # end for: i
v_ints <- unlist(l_ints)
v_ints <- v_ints[!(v_ints==0)]
v_ints2 <- v_ints[v_ints<11] # subset

# Plot
par(mfrow=c(1,2))
hist(v_ints, breaks=seq(0, 50, length=50), xlab="Interval Size", main="")
barplot(table(v_ints2), xlab="Interval Size", ylab="Frequency")


# ----------------------------------------------------------
# -------- Looking at Data ---------------------------------
# ----------------------------------------------------------

# --------- Line graphs & Marginals --------

# Subset 1 person
i <- 1
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id==u_subj[i], ]
X <- data_i$happy
# Layout
lmat <- matrix(1:2, nrow=1)
layout(lmat, widths = c(1, .3))
# Time Series
par(mar=c(4,4,2,1))
plot.new()
plot.window(xlim=c(1,240), ylim=c(0,100))
axis(1, c(0,50,100, 150, 200, 240))
axis(2, las=2)
abline(v=seq(0, 240, length=41), col="lightgrey", lty=2) # day indicator
lines(X, col="black")
title(xlab="Time", ylab="Response", main="Variable 'Happy'; Subject 1", font.main=1)
# Marginals
par(mar=c(4,0,2,.5))
breaks <- seq(0, 100, length=30)
tb <- hist(X, breaks=breaks, plot = FALSE)
barplot(tb$counts, axes = FALSE, horiz = TRUE, ylab="")

# --------- Line graphs & Marginals: All individuals --------

# Same as above, but for all data, with 1 page=1subject

pdf("figures/Looking1_tsandmarginals.pdf", width=8, height=8)

u_subj <- unique(data$subj_id)
names <- colnames(data)[5:12]

# Layout
lmat <- matrix(1:16, nrow=4, byrow = TRUE)
lo <- layout(lmat, widths = c(1,.3,1,.3))

# Loop through data
for(i in 1:125) {

  data_i <- data[data$subj_id==u_subj[i], ]

  for(j in 1:8) {

    X <- data_i[, 4+j]
    # Time Series
    par(mar=c(4,4,2,1))
    plot.new()
    plot.window(xlim=c(1,240), ylim=c(0,100))
    axis(1, c(0,50,100, 150, 200, 240))
    axis(2, las=2)
    abline(v=seq(0, 240, length=41), col="lightgrey", lty=2) # day indicator
    lines(X, col="black")
    title(xlab="Time", ylab="Response", line=2)
    title(main=paste0("item: ", names[j], "; subj: ", i), font.main=1)
    # Marginals
    par(mar=c(4,0,2,.5))
    breaks <- seq(0, 100, length=30)
    tb <- hist(X, breaks=breaks, plot = FALSE)
    barplot(tb$counts, axes = FALSE, horiz = TRUE, ylab="")

  } # end for: j


} # end for: i


dev.off()


# --------- Bivariate Plot: Happy & Sad --------

library(scales)

# Subset 1 person


i <- 2
u_subj <- unique(data$subj_id)
par(mar=c(4,4,1,1))
plot.new()
plot.window(xlim=c(0,100), ylim=c(0,100))
axis(1)
axis(2, las=2)
title(xlab="Happy", ylab="Sad", line=2.5, cex.lab=1.5)
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id==u_subj[i], ]
points(data_i$happy, data_i$sad, pch=20, col=alpha("black", alpha=.5), cex=2)
title(main=paste0("Subject: ", i), font.main = 1)


# --------- Do the same for each subject --------

pdf("figures/Looking2_biv_happy_sad.pdf", width=9, height=9)

par(mfrow=c(3,3))

u_subj <- unique(data$subj_id)
par(mar=c(4,4,2,1))

for(i in 1:125){

  # Layout
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0,100))
  axis(1)
  axis(2, las=2)
  title(xlab="Happy", ylab="Sad", line=2.5)
  title(main=paste0("Subject: ", i), font.main = 1)

  # Data
  data_i <- data[data$subj_id==u_subj[i], ]
  data_i <- data[data$subj_id==u_subj[i], ]
  points(data_i$happy, data_i$sad, pch=20, col=alpha("black", alpha=.5), cex=1.5)

}

dev.off()



# ----------------------------------------------------------
# -------- Descriptives ------------------------------------
# ----------------------------------------------------------

# -------- Mean & SD ---------

# Subset 1 person
i <- 1
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id==u_subj[i], ]
X <- data_i$happy
# Layout
lmat <- matrix(1:2, nrow=1)
layout(lmat, widths = c(1, .3))
# Time Series
par(mar=c(4,4,2,1))
plot.new()
plot.window(xlim=c(1,240), ylim=c(0,100))
axis(1, c(0,50,100, 150, 200, 240))
axis(2, las=2)
abline(v=seq(0, 240, length=41), col="lightgrey", lty=2) # day indicator
lines(X, col="black")
title(xlab="Time", ylab="Response", main="Variable 'Happy'; Subject 1", font.main=1)
# Marginals
par(mar=c(4,0,2,.5))
breaks <- seq(0, 100, length=30)
tb <- hist(X, breaks=breaks, plot = FALSE)
barplot(tb$counts, axes = FALSE, horiz = TRUE, ylab="")
sc <-  (100/29) #
X_sc <- X / sc
mean_X <- mean(X_sc, na.rm=TRUE) # scale to barplot
sd_X <- sd(X_sc, na.rm=TRUE)
abline(h=mean_X, col="blue", lwd=2)
segments(x0 = 15, y0 = mean_X-5, x1=15, y1 = mean_X+5, lwd=2, lty=3, col="blue")
segments(x0 = 10, y0 = mean_X-5, x1=20, y1 = mean_X-5, lwd=2, lty=3, col="blue")
segments(x0 = 10, y0 = mean_X+5, x1=20, y1 = mean_X+5, lwd=2, lty=3, col="blue")


# -------- Mean & SD: Heterogeneity ---------

# --- Compute means/sds ---
m_res <- matrix(NA, 125, 4) # means (sad,happy), sds (sad, happy)
u_subj <- unique(data$subj_id)
for(i in 1:125) {
  data_i <- data[data$subj_id==u_subj[i], ]
  m_res[i, 1] <- mean(data_i$happy, na.rm=TRUE)
  m_res[i, 2] <- mean(data_i$sad, na.rm=TRUE)
  m_res[i, 3] <- sd(data_i$happy, na.rm=TRUE)
  m_res[i, 4] <- sd(data_i$sad, na.rm=TRUE)
}
# --- Figure ---
lmat <- matrix(5:8, nrow=2, byrow = TRUE)
lmat <- cbind(1:2, lmat)
lmat <- rbind(c(0,3,4), lmat)
lo <- layout(lmat, widths = c(.15, 1, 1), heights = c(0.15, 1, 1))
# Plot Labels
PL("Mean", srt = 90, cex=1.5)
PL("SD", srt = 90, cex=1.5)
PL("Happy", cex=1.5)
PL("Sad", cex=1.5)
# Plot Data
par(mar=c(3,3,2,1))
breaks <- seq(0, 100, length=30)
for(i in 1:4) hist(m_res[, i], main="", xlim=c(0,100),
                   xlab="", ylab="", breaks=breaks)

# -------- Skewness & Modes ---------

par(mfrow=c(1,2))
i <- 1
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id==u_subj[i], ]
breaks <- seq(0, 100, length=30)
hist(data_i$happy, xlim=c(0,100),
     xlab="", ylab="", breaks=breaks, main="Happy (Subj 1)", font.main=1)
i <- 5
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id==u_subj[i], ]
hist(data_i$anxious, xlim=c(0,100),
     xlab="", ylab="", breaks=breaks, main="Anxious (Subj 5)", font.main=1)


# -------- RMSE & AR ---------

## Make four examples
# Functions for AR/RMSSD
AR <- function(x) {
  cor(x[-1], x[-length(x)])
}
RMSSD <- function(x) {
  sqrt(mean((x[-1] -x[-length(x)])^2))
}
set.seed(6)
x <- 1:240
# Example 1: step function
y1 <- c(rep(40, 120), rep(60, 120))
# Example 2: Gaussian noise [matched]
y2 <- rnorm(240, 50, 10)
# Example 3: Gaussian noise [matched]
y3 <- rep(40, 240)
for(i in 2:240) y3[i] <- 20 + 0.65*y3[i-1] + rnorm(1, 0, 5)
# Example 2: Gaussian noise [matched]
y4 <- rep(50, 240)
for(i in 2:240) y4[i] <- 1*y4[i-1] + rnorm(1, 0, 3)
l_y <- list(y1, y2, y3, y4)
# Plotting
par(mfrow=c(2,2), mar=c(4,3,2,1))
for(i in 1:4) {
  plot.new()
  plot.window(xlim=c(1,240), ylim=c(0,100))
  axis(1, c(0,50,100, 150, 200, 240))
  axis(2, las=2)
  cex <- 1.05
  text(175, 100, paste0("SD = ", round(sd(l_y[[i]]),2)), adj=0, cex=cex)
  text(175, 90, paste0("AR = ", round(AR(l_y[[i]]),2)), adj=0, cex=cex)
  text(175, 80, paste0("RMSSD = ", round(RMSSD(l_y[[i]]),2)), adj=0, cex=cex)
  lines(x, l_y[[i]], lwd=2)
}

# ----------------------------------------------------------
# -------- Basic Modeling: AR ------------------------------
# ----------------------------------------------------------

# ----- Generate Data from AR model -----

out <- lm(X[-1] ~ X[-length(X)])
out
# Generate Data from AR model
set.seed(1)
X_sim <- rep(NA, 240)
X_sim[1] <- X[1]
coefs <- out$coefficients
for(i in 2:240) X_sim[i] <- coefs[1] + coefs[2]*X_sim[i-1] + rnorm(1, 0, sd(residuals(out)))

# Subset 1 person
i <- 1
u_subj <- unique(data$subj_id)
data_i <- data[data$subj_id==u_subj[i], ]
X <- data_i$happy
# Layout
lmat <- matrix(1:4, nrow=1)
layout(lmat, widths = c(1, .3, 1, .3))

TSplotARfit <- function(X, main=NULL) {
  # Time Series
  par(mar=c(4,4,2,1))
  plot.new()
  plot.window(xlim=c(1,240), ylim=c(-10,130))
  axis(1, c(0,50,100, 150, 200, 240))
  axis(2, las=2)
  abline(v=seq(0, 240, length=41), col="lightgrey", lty=2) # day indicator
  lines(X, col="black")
  title(xlab="Time", ylab="Response", main=main, font.main=1)
  # Marginals
  par(mar=c(4,0,2,.5))
  breaks <- seq(-10, 130, length=30)
  tb <- hist(X, breaks=breaks, plot = FALSE)
  barplot(tb$counts, axes = FALSE, horiz = TRUE, ylab="")
}
TSplotARfit(X, main="Variable 'Happy'; Subject 1")
TSplotARfit(X_sim, main="Simulated data")



# ----- Fit VAR model to Subj1 -----

library(mgm)
out <- mvar(data=data_i_noNA_vars,
            type=rep("g", 8),
            level=rep(1, 8),
            lags = 1,
            lambdaSel = "EBIC",
            lambdaSeq = 0,
            dayvar = data_i_noNA$dayno,
            beepvar = data_i_noNA$beep,
            threshold = "none",
            scale=FALSE)

# Phi-matrix
out$signs[is.na(out$signs)] <- 1
phi <- out$wadj[, , 1] * out$signs[, , 1]




# ------------------------------------------------------
# -------- Basic Modeling: VAR -------------------------
# ------------------------------------------------------

# ----- Select 2 PE, 2 NE -----
data_i <- data[data$subj_id==1, ]
data_i_vars <- data_i[, c("happy", "relaxed", "anxious", "sad")]
data_i_noNA <- na.omit(data_i)
data_i_noNA_vars <- data_i_noNA[, c("happy", "relaxed", "anxious", "sad")]

# ----- Show data -----
# Layout
lmat <- matrix(1:8, nrow=2, byrow = TRUE)
lo <- layout(lmat, widths = c(1, .3, 1, .3))
# layout.show(lo)
TSplotARfit(data_i$happy, main="Happy")
TSplotARfit(data_i$relaxed, main="Relaxed")
TSplotARfit(data_i$anxious, main="anxious")
TSplotARfit(data_i$sad, main="Sad")


# ----- Fit VAR -----
out <- mvar(data=data_i_noNA_vars,
            type=rep("g", 4),
            level=rep(1, 4),
            lags = 1,
            lambdaSel = "EBIC",
            lambdaSeq = 0,
            dayvar = data_i_noNA$dayno,
            beepvar = data_i_noNA$beep,
            threshold = "none",
            scale=FALSE)

# intercepts
ints <- unlist(out$intercepts)

# phi-matrix
out$signs[is.na(out$signs)] <- 1
phi <- out$wadj[, , 1] * out$signs[, , 1]

# residual variances
pred <- predict(out, data=data_i_noNA_vars,
                dayvar = data_i_noNA$dayno,
                beepvar = data_i_noNA$beep)
res <- pred$predicted - data_i_noNA_vars
cov_res <- cov(res)

# ----- Make Graph -----

qgraph(phi, labels=colnames(data_i_noNA_vars), edge.labels=T)


# ----- Generate Data -----

set.seed(14)
data_sim <- simulateVAR(pars = phi, means = ints,
                        Nt = 240,
                        init = as.numeric(data_i_noNA_vars[1,]),
                        residuals = cov_res) # inital values
range(unlist(data_sim))

# Combine into array
l_data <- list(data_i_vars, data_sim)

# ----- Plot Figure -----
par(mfrow=c(1,2))

# Empirical data
lim <- c(-20,120)
par(mar=c(4,4,1,1))
plot.new()
plot.window(xlim=lim, ylim=lim)
axis(1)
axis(2, las=2)
abline(h=c(0,100), col="grey")
abline(v=c(0,100), col="grey")
title(xlab="Happy", ylab="Sad", line=2.5, cex.lab=1.5)
points(data_i$happy, data_i$sad, pch=20, col=alpha("black", alpha=.5), cex=2)
title(main=paste0("Subject: 1"), font.main = 1)

# Simulated data
par(mar=c(4,4,1,1))
plot.new()
plot.window(xlim=lim, ylim=lim)
axis(1)
axis(2, las=2)
abline(h=c(0,100), col="grey")
abline(v=c(0,100), col="grey")
title(xlab="Happy", ylab="Sad", line=2.5, cex.lab=1.5)
points(data_sim[, 1], data_sim[, 2], pch=20, col=alpha("black", alpha=.5), cex=2)
title(main="Simulated Data", font.main = 1)











