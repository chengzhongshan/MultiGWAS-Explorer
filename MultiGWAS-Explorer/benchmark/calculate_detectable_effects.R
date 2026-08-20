args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript calculate_detectable_effects.R input.tsv output.tsv")
d <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)

alphas <- c(`genome_wide_5e-8` = 5e-8, `exploratory_1e-5` = 1e-5)
target_power <- 0.80
z_power <- qnorm(target_power)

power_two_sided <- function(beta, se, alpha) {
  zcrit <- qnorm(1 - alpha / 2)
  ncp <- abs(beta) / se
  pnorm(-zcrit - ncp) + 1 - pnorm(zcrit - ncp)
}

rows <- list()
k <- 0L
for (i in seq_len(nrow(d))) {
  maf_f <- min(d$EAF_FEMALE[i], 1 - d$EAF_FEMALE[i])
  maf_m <- min(d$EAF_MALE[i], 1 - d$EAF_MALE[i])
  beta_diff <- d$BETA_FEMALE[i] - d$BETA_MALE[i]
  se_diff <- sqrt(d$SE_FEMALE[i]^2 + d$SE_MALE[i]^2)
  tests <- data.frame(
    TEST = c("FEMALE_LOG_OR", "MALE_LOG_OR", "FEMALE_MINUS_MALE_LOG_OR"),
    OBSERVED_BETA = c(d$BETA_FEMALE[i], d$BETA_MALE[i], beta_diff),
    SE = c(d$SE_FEMALE[i], d$SE_MALE[i], se_diff),
    MAF = c(maf_f, maf_m, min(maf_f, maf_m)),
    stringsAsFactors = FALSE
  )
  for (aname in names(alphas)) {
    alpha <- alphas[[aname]]
    zcrit <- qnorm(1 - alpha / 2)
    for (j in seq_len(nrow(tests))) {
      k <- k + 1L
      mde <- (zcrit + z_power) * tests$SE[j]
      rows[[k]] <- data.frame(
        SNP = d$SNP[i], PAIR = d$PAIR[i], A1 = d$A1[i], A2 = d$A2[i],
        TEST = tests$TEST[j], ALPHA_LABEL = aname, ALPHA = alpha,
        TARGET_POWER = target_power, MAF = tests$MAF[j],
        OBSERVED_BETA = tests$OBSERVED_BETA[j], SE = tests$SE[j],
        OBSERVED_OR_OR_RATIO = exp(tests$OBSERVED_BETA[j]),
        MIN_DETECTABLE_ABS_BETA = mde,
        MIN_DETECTABLE_OR_OR_RATIO_LOWER = exp(-mde),
        MIN_DETECTABLE_OR_OR_RATIO = exp(mde),
        POWER_AT_OBSERVED_EFFECT = power_two_sided(tests$OBSERVED_BETA[j], tests$SE[j], alpha),
        NCASE_FEMALE = d$NCASE_FEMALE[i], NCONTROL_FEMALE = d$NCONTROL_FEMALE[i],
        NCASE_MALE = d$NCASE_MALE[i], NCONTROL_MALE = d$NCONTROL_MALE[i],
        INFO_FEMALE = d$INFO_FEMALE[i], INFO_MALE = d$INFO_MALE[i],
        ASSUMPTIONS = "normal-Wald; independent strata; observed-SE-based; beta is log(OR)",
        stringsAsFactors = FALSE
      )
    }
  }
}
out <- do.call(rbind, rows)
write.table(out, args[2], sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
