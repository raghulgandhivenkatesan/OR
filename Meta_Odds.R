#################### Load necessary libraries #######################
library(metafor)
library(meta)

#################### Data for Odds Ratio Meta-Analysis ##############
# Replace with your actual data

or_data <- data.frame(
  Author = c("Smith 2020", "Lee 2021", "Chen 2022"),
  ai = c(10, 15, 5), # Events in treatment (dengue with myocarditis)
  n1i = c(100, 120, 90),   # Total in treatment group
  ci = c(2, 5, 3),    # Events in control (dengue without myocarditis)
  n2i = c(100, 120, 90)     # Total in control group
)

#################### Step 1: Calculate log odds ratios ######################
or_escalc <- escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = or_data)
print(or_escalc)

#################### Step 2: Random-effects meta-analysis ##################
or_rma <- rma(yi, vi, data = or_escalc, method = "DL", weighted = TRUE)
print(or_rma, digits = 2)
confint(or_rma)

#################### Step 3: Identify outliers via residuals ################
or_resid <- rstudent(or_rma)
abs_z_or <- abs(or_resid$z)
outliers_or <- or_resid[order(-abs_z_or), ]
print(outliers_or)

#################### Step 4: Leave-One-Out Analysis ########################
or_l1o <- leave1out(or_rma)
yi_l1o <- or_l1o$estimate
sei_l1o <- or_l1o$se
ci.lb <- yi_l1o - 1.96 * sei_l1o
ci.ub <- yi_l1o + 1.96 * sei_l1o

forest(yi_l1o, sei = sei_l1o,
       slab = or_data$Author,
       xlab = "Log Odds Ratio (Leave-One-Out)",
       refline = or_rma$b,
       digits = 3,
       alim = c(min(ci.lb), max(ci.ub)),
       cex = 0.8,
       col = "blue",
       main = "Leave-One-Out Analysis (OR)")

#################### Step 5: Baujat Plot ##################################
baujat(or_rma, main = "Baujat Plot for Heterogeneity (OR)")

#################### Step 6: Influence Diagnostics #########################
or_inf <- influence(or_rma)
print(or_inf)
plot(or_inf)

#################### Step 7: Remove Outliers (if needed) ##################
or_no_out <- or_escalc[-c(1, 6), ]  # Replace with detected outlier indices
or_rma_updated <- rma(yi, vi, data = or_no_out, method = "DL", weighted = TRUE)
print(or_rma_updated, digits = 2)

#################### Step 8: Forest Plot (OR) #############################
meta_or <- metabin(
  event.e = ai,
  n.e = n1i,
  event.c = ci,
  n.c = n2i,
  studlab = Author,
  data = or_data,
  sm = "OR", method = "Inverse", method.tau = "DL", 
  incr = 0.5, allstudies = TRUE, hakn = TRUE
)

png("forestplot_or.png", width = 1000, height = 1000)
forest(meta_or,
       xlim = c(0.01, 100), # Adjust based on your data
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("OR", "95% C.I.", "Weight"),
       leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c"),
       leftlabs = c("Study", "Events (Dengue+Myocarditis)", "N", "Events (Control)", "N"),
       squaresize = 0.5,
       col.square = "darkgreen",
       col.diamond = "red",
       col.diamond.lines = "black",
       print.I2 = TRUE, print.Q = TRUE,
       main = "Odds Ratio Meta-Analysis")
dev.off()


#################### Step 9: Funnel Plot ##################################
funnel(or_rma, main = "Funnel Plot (Odds Ratio)")

#################### Step 10: Trim-and-Fill for Bias ######################
or_trimfill <- trimfill(or_rma)
print(or_trimfill)
funnel(or_trimfill, main = "Trim-and-Fill Funnel Plot (OR)")

#################### Step 11: Egger's Test ################################
egger_or <- regtest(or_rma, model = "lm", predictor = "sei")
print(egger_or)

#################### Step 12: Rank Correlation Test #######################
rank_corr_or <- ranktest(or_rma)
print(rank_corr_or)

#################### Step 13: Meta-regression by Sample Size #############
# Total sample per study = n1i + n2i
or_escalc$N_total <- or_data$n1i + or_data$n2i
meta_reg_or <- rma(yi, vi, mods = ~ N_total, data = or_escalc, method = "DL")
summary(meta_reg_or)

#################### Step 14: Bubble Plot for Meta-Regression ############
bubble_or <- predict(meta_reg_or, newmods = or_escalc$N_total)
plot(or_escalc$N_total, or_escalc$yi,
     xlab = "Sample Size", 
     ylab = "Log Odds Ratio", 
     main = "Bubble Plot (Meta-Regression by Sample Size)",
     pch = 21, bg = "lightgray", cex = 1.2)
lines(or_escalc$N_total, bubble_or$pred, col = "red", lwd = 2)
legend("topright", legend = "Regression Line", col = "red", lwd = 2)

#################### Step 15: Distribution of True Effect ################
tau2_plot <- ifelse(or_rma$tau2 == 0, 0.01, or_rma$tau2)

true_or <- seq(or_rma$b - 4 * sqrt(tau2_plot), or_rma$b + 4 * sqrt(tau2_plot), length.out = 1000)
density_or <- dnorm(true_or, mean = or_rma$b, sd = sqrt(tau2_plot))

plot(true_or, density_or, type = "l", lwd = 2, col = "blue",
     xlab = "True Log Odds Ratio", ylab = "Density", main = "Distribution of the True Effect")
abline(v = or_rma$b, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Distribution", "Overall Log OR"), col = c("blue", "red"), lty = c(1, 2), lwd = 2)