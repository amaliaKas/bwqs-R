##############################################
##########Project: HELIX study - adolescence - PFAS mixtures and liver outcomes
##########Creator: Amalia Kasmiridou
##########Aim: BWQS main + sensitivity Analyses for imp=1
##########Last update: 22/07/2026 
##############################################
# Set working directory & load the dataset ####
setwd("~/BWQS/")
load("data_pfas_liver.Rdata")

# Imputed dataset #1
data1 <- subset(data_pfas_liver, .imp==1)

# Datasets for sensitivity analysis among participants with prenatal PFAS 
# Define the variables we want to check for NAs
# _m refers to prenatal PFAS, _c refers to childhood PFAS, _a refers to adolescent PFAS
varsc <- c("log2_pfoa_m", "log2_pfna_m", "log2_pfhxs_m","log2_pfos_m",
          "log2_pfoa_c", "log2_pfna_c", "log2_pfhxs_c","log2_pfos_c", 
          "log2_pfunda_c")
data1c <- data1[complete.cases(data1[, varsc]), ]
varsa <- c("log2_pfoa_m", "log2_pfna_m", "log2_pfhxs_m","log2_pfos_m",
           "log2_pfoa_a", "log2_pfna_a", "log2_pfhxs_a","log2_pfos_a",
           "log2_pfunda_a", "log2_pfda_a", "log2_pfhpa_a", "log2_x9clpfesa_a")
data1a <- data1[complete.cases(data1[, varsa]), ]

# BWQS package ####
install.packages("remotes")
remotes::install_github("ElenaColicino/bwqs")
library(BWQS)

# Define Mixtures ####
premix <- c("log2_pfoa_m", "log2_pfna_m", "log2_pfhxs_m", "log2_pfos_m")
childmix <- c("log2_pfoa_c", "log2_pfna_c", "log2_pfhxs_c", "log2_pfos_c", 
              "log2_pfunda_c")
adolmix1 <- c("log2_pfoa_a", "log2_pfna_a", "log2_pfhxs_a","log2_pfos_a")
adolmix2 <- c("log2_pfoa_a", "log2_pfna_a", "log2_pfhxs_a","log2_pfos_a",
              "log2_pfunda_a", "log2_pfda_a", "log2_pfhpa_a", "log2_x9clpfesa_a")
##############################################
# MAIN MODEL ####
# 1. ALT ####

# Prenatal
bwqs.alt.1m <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                      native + fishpreg_tert,
                    mix_name = premix,
                    data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output 
bwqs.alt.1m_fit <- bwqs.alt.1m$summary_fit
save(bwqs.alt.1m_fit, file = "main/bwqs.alt.1m.Rdata")
write.xlsx(bwqs.alt.1m_fit, "main/bwqs.alt.1m.xlsx")

# Childhood
bwqs.alt.1c <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                      bw + bf + fishchild_tert + ethn2c,
                    mix_name = childmix,
                    data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
bwqs.alt.1c_fit <- bwqs.alt.1c$summary_fit
## model output 
bwqs.alt.1c_fit <- bwqs.alt.1c$summary_fit
save(bwqs.alt.1c_fit, file = "main/bwqs.alt.1c.Rdata")
write.xlsx(bwqs.alt.1c_fit, "main/bwqs.alt.1c.xlsx")

# Adolescence Mix 1
bwqs.alt.1a.1 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                        bw + bf + fishadol_tert + ethn2c,
                      mix_name = adolmix1,
                      data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.1a.1_fit <- bwqs.alt.1a.1$summary_fit
save(bwqs.alt.1a.1_fit, file = "main/bwqs.alt.1a.1.Rdata")
write.xlsx(bwqs.alt.1a.1_fit, "main/bwqs.alt.1a.1.xlsx")

# Adolescence Mix 2
bwqs.alt.1a.2 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                        bw + bf + fishadol_tert + ethn2c,
                      mix_name = adolmix2,
                      data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.1a.2_fit <- bwqs.alt.1a.2$summary_fit
save(bwqs.alt.1a.2_fit, file = "main/bwqs.alt.1a.2.Rdata")
write.xlsx(bwqs.alt.1a.2_fit, "main/bwqs.alt.1a.2.xlsx")

# 2. MASLD ####

# Prenatal
bwqs.masld.1m <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                        native + fishpreg_tert,
                      mix_name = premix,
                      data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output 
bwqs.masld.1m_fit <- bwqs.masld.1m$summary_fit
save(bwqs.masld.1m_fit, file = "main/bwqs.masld.1m.Rdata")
write.xlsx(bwqs.masld.1m_fit, "main/bwqs.masld.1m.xlsx")

# Childhood
bwqs.masld.1c <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                        bw + bf + fishchild_tert + ethn2c,
                      mix_name = childmix
                      data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
bwqs.masld.1c_fit <- bwqs.masld.1c$summary_fit
## model output 
bwqs.masld.1c_fit <- bwqs.masld.1c$summary_fit
save(bwqs.masld.1c_fit, file = "main/bwqs.masld.1c.Rdata")
write.xlsx(bwqs.masld.1c_fit, "main/bwqs.masld.1c.xlsx")

# Adolescence Mix 1
bwqs.masld.1a.1 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                          bw + bf + fishadol_tert + ethn2c,
                        mix_name = adolmix1,
                        data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.1a.1_fit <- bwqs.masld.1a.1$summary_fit
save(bwqs.masld.1a.1_fit, file = "main/bwqs.masld.1a.1.Rdata")
write.xlsx(bwqs.masld.1a.1_fit, "main/bwqs.masld.1a.1.xlsx")

# Adolescence Mix 2
bwqs.masld.1a.2 <- bwqs(masld_2c ~ cohort + sex + mbmi + age + edumc + parity + 
                          globalsmok_m + visit_age_years + weights_atr +
                          bw + bf + fishadol_tert + ethn2c,
                        mix_name = adolmix2,
                        data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.1a.2_fit <- bwqs.masld.1a.2$summary_fit
save(bwqs.masld.1a.2_fit, file = "main/bwqs.masld.1a.2.Rdata")
write.xlsx(bwqs.masld.1a.2_fit, "main/bwqs.masld.1a.2.xlsx")

##############################################
# SENSITIVITY MODEL #1 (population with prenatal PFAS) ####
# 1. ALT ####

# Childhood
bwqs.alt.sens.1c <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                           bw + bf + fishchild_tert + ethn2c,
                         mix_name = childmix,
                         data = data1c, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output 
bwqs.alt.sens.1c_fit <- bwqs.alt.sens.1c$summary_fit
save(bwqs.alt.sens.1c_fit, file = "sens/bwqs.alt.sens.1c.Rdata")
write.xlsx(bwqs.alt.sens.1c_fit, "sens/bwqs.alt.sens.1c.xlsx")

# Adolescence Mix 1
bwqs.alt.sens.1a.1 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             bw + bf + fishadol_tert + ethn2c,
                           mix_name = adolmix1,
                           data = data1a, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.sens.1a.1_fit <- bwqs.alt.sens.1a.1$summary_fit
save(bwqs.alt.sens.1a.1_fit, file = "sens/bwqs.alt.sens.1a.1.Rdata")
write.xlsx(bwqs.alt.sens.1a.1_fit, "sens/bwqs.alt.sens.1a.1.xlsx")

# Adolescence Mix 2
bwqs.alt.sens.1a.2 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             bw + bf + fishadol_tert + ethn2c,
                           mix_name = adolmix2,
                           data = data1a, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.sens.1a.2_fit <- bwqs.alt.sens.1a.2$summary_fit
save(bwqs.alt.sens.1a.2_fit, file = "sens/bwqs.alt.sens.1a.2.Rdata")
write.xlsx(bwqs.alt.sens.1a.2_fit, "sens/bwqs.alt.sens.1a.2.xlsx")

# 2. MASLD ####

# Childhood
bwqs.masld.sens.1c <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             bw + bf + fishchild_tert + ethn2c,
                           mix_name = childmix,
                           data = data1c, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
bwqs.masld.sens.1c_fit <- bwqs.masld.sens.1c$summary_fit
## model output 
bwqs.masld.sens.1c_fit <- bwqs.masld.sens.1c$summary_fit
save(bwqs.masld.sens.1c_fit, file = "sens/bwqs.masld.sens.1c.Rdata")
write.xlsx(bwqs.masld.sens.1c_fit, "sens/bwqs.masld.sens.1c.xlsx")

# Adolescence Mix 1
bwqs.masld.sens.1a.1 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                               bw + bf + fishadol_tert + ethn2c,
                             mix_name = adolmix1,
                             data = data1a, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.sens.1a.1_fit <- bwqs.masld.sens.1a.1$summary_fit
save(bwqs.masld.sens.1a.1_fit, file = "sens/bwqs.masld.sens.1a.1.Rdata")
write.xlsx(bwqs.masld.sens.1a.1_fit, "sens/bwqs.masld.sens.1a.1.xlsx")

# Adolescence Mix 2
bwqs.masld.sens.1a.2 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                               bw + bf + fishadol_tert + ethn2c,
                             mix_name = adolmix2,
                             data = data1a, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.sens.1a.2_fit <- bwqs.masld.sens.1a.2$summary_fit
save(bwqs.masld.sens.1a.2_fit, file = "sens/bwqs.masld.sens.1a.2.Rdata")
write.xlsx(bwqs.masld.sens.1a.2_fit, "sens/bwqs.masld.sens.1a.2.xlsx")

##############################################
# SENSITIVITY MODEL #2 (Adjustment for previous PFAS exposures) ####
# 1. ALT ####

# Childhood
bwqs.alt.1c_m <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                        bw + bf + fishchild_tert + ethn2c + 
                        log2_pfoa_m + log2_pfna_m + log2_pfhxs_m + log2_pfos_m,
                      mix_name = childmix,
                      data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output 
bwqs.alt.1c_m_fit <- bwqs.alt.1c_m$summary_fit
save(bwqs.alt.1c_m_fit, file = "prevadj/bwqs.alt.1c_m.Rdata")
write.xlsx(bwqs.alt.1c_m_fit, "prevadj/bwqs.alt.1c_m.xlsx")

# Adolescence Mix 1
bwqs.alt.1a.1_m <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                          bw + bf + fishadol_tert + ethn2c + 
                          log2_pfoa_m + log2_pfna_m + log2_pfhxs_m + log2_pfos_m + 
                          log2_pfoa_c + log2_pfna_c + log2_pfhxs_c + log2_pfos_c + log2_pfunda_c,
                        mix_name = adolmix1,
                        data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.1a.1_m_fit <- bwqs.alt.1a.1_m$summary_fit
save(bwqs.alt.1a.1_m_fit, file = "prevadj/bwqs.alt.1a.1_m-2.Rdata")
write.xlsx(bwqs.alt.1a.1_m_fit, "prevadj/bwqs.alt.1a.1_m-2.xlsx")

# Adolescence Mix 2
bwqs.alt.1a.2_m <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                          bw + bf + fishadol_tert + ethn2c + 
                          log2_pfoa_m + log2_pfna_m + log2_pfhxs_m + log2_pfos_m + 
                          log2_pfoa_c + log2_pfna_c + log2_pfhxs_c + log2_pfos_c + log2_pfunda_c,
                        mix_name = adolmix2,
                        data = data1, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.1a.2_m_fit <- bwqs.alt.1a.2_m$summary_fit
save(bwqs.alt.1a.2_m_fit, file = "prevadj/bwqs.alt.1a.2_m-2.Rdata")
write.xlsx(bwqs.alt.1a.2_m_fit, "prevadj/bwqs.alt.1a.2_m-2.xlsx")

# 2. MASLD ####

# Childhood
bwqs.masld.1c_m <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                          bw + bf + fishchild_tert + ethn2c+ 
                          log2_pfoa_m + log2_pfna_m + log2_pfhxs_m + log2_pfos_m,
                        mix_name = childmix,
                        data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output 
bwqs.masld.1c_m_fit <- bwqs.masld.1c_m$summary_fit
save(bwqs.masld.1c_m_fit, file = "prevadj/bwqs.masld.1c_m.Rdata")
write.xlsx(bwqs.masld.1c_m_fit, "prevadj/bwqs.masld.1c_m.xlsx")

# Adolescence Mix 1
bwqs.masld.1a.1_m <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                            bw + bf + fishadol_tert + ethn2c + 
                            log2_pfoa_m + log2_pfna_m + log2_pfhxs_m + log2_pfos_m + 
                            log2_pfoa_c + log2_pfna_c + log2_pfhxs_c + log2_pfos_c + log2_pfunda_c,
                          mix_name = adolmix1,
                          data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.1a.1_m_fit <- bwqs.masld.1a.1_m$summary_fit
save(bwqs.masld.1a.1_m_fit, file = "prevadj/bwqs.masld.1a.1_m-2.Rdata")
write.xlsx(bwqs.masld.1a.1_m_fit, "prevadj/bwqs.masld.1a.1_m-2.xlsx")

# Adolescence Mix 2
bwqs.masld.1a.2_m <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                            bw + bf + fishadol_tert + ethn2c + 
                            log2_pfoa_m + log2_pfna_m + log2_pfhxs_m + log2_pfos_m + 
                            log2_pfoa_c + log2_pfna_c + log2_pfhxs_c + log2_pfos_c + log2_pfunda_c,
                          mix_name = adolmix2,
                          data = data1, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.1a.2_m_fit <- bwqs.masld.1a.2_m$summary_fit
save(bwqs.masld.1a.2_m_fit, file = "prevadj/bwqs.masld.1a.2_m-2.Rdata")
write.xlsx(bwqs.masld.1a.2_m_fit, "prevadj/bwqs.masld.1a.2_m-2.xlsx")

##############################################
# SENSITIVITY MODEL #3 (population w/ PGS data) ####
data1.pgs <- subset(data1, !is.na(pgs2071_2c)|!is.na(pgs0027_2c))

# 1. ALT ####

# Prenatal
bwqs.alt.pgs.1m <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                          native + fishpreg_tert,
                        mix_name = premix,
                        data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output 
bwqs.alt.pgs.1m_fit <- bwqs.alt.pgs.1m$summary_fit
save(bwqs.alt.pgs.1m_fit, file = "pgs/bwqs.alt.pgs.1m.Rdata")
write.xlsx(bwqs.alt.pgs.1m_fit, "pgs/bwqs.alt.pgs.1m.xlsx")

# Childhood
bwqs.alt.pgs.1c <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                          bw + bf + fishchild_tert + ethn2c,
                        mix_name = childmix,
                        data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
bwqs.alt.pgs.1c_fit <- bwqs.alt.pgs.1c$summary_fit
## model output 
bwqs.alt.pgs.1c_fit <- bwqs.alt.pgs.1c$summary_fit
save(bwqs.alt.pgs.1c_fit, file = "pgs/bwqs.alt.pgs.1c.Rdata")
write.xlsx(bwqs.alt.pgs.1c_fit, "pgs/bwqs.alt.pgs.1c.xlsx")

# Adolescence Mix 1
bwqs.alt.pgs.1a.1 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                            bw + bf + fishadol_tert + ethn2c,
                          mix_name = adolmix1,
                          data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.pgs.1a.1_fit <- bwqs.alt.pgs.1a.1$summary_fit
save(bwqs.alt.pgs.1a.1_fit, file = "pgs/bwqs.alt.pgs.1a.1.Rdata")
write.xlsx(bwqs.alt.pgs.1a.1_fit, "pgs/bwqs.alt.pgs.1a.1.xlsx")

# Adolescence Mix 2
bwqs.alt.pgs.1a.2 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                            bw + bf + fishadol_tert + ethn2c,
                          mix_name = adolmix2,
                          data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.alt.pgs.1a.2_fit <- bwqs.alt.pgs.1a.2$summary_fit
save(bwqs.alt.pgs.1a.2_fit, file = "pgs/bwqs.alt.pgs.1a.2.Rdata")
write.xlsx(bwqs.alt.pgs.1a.2_fit, "pgs/bwqs.alt.pgs.1a.2.xlsx")

# 2. MASLD ####

# Prenatal
bwqs.masld.pgs.1m <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                            native + fishpreg_tert,
                          mix_name = premix,
                          data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
bwqs.masld.pgs.1m_fit <- bwqs.masld.pgs.1m$summary_fit
## model output 
bwqs.masld.pgs.1m_fit <- bwqs.masld.pgs.1m$summary_fit
save(bwqs.masld.pgs.1m_fit, file = "P:/HELIX chemicals/PFAS/LIVER/BWQS/imp1/pgs/bwqs.masld.pgs.1m.Rdata")
write.xlsx(bwqs.masld.pgs.1m_fit, "P:/HELIX chemicals/PFAS/LIVER/BWQS/imp1/pgs/bwqs.masld.pgs.1m.xlsx")


# Childhood
bwqs.masld.pgs.1c <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                            bw + bf + fishchild_tert + ethn2c,
                          mix_name = childmix,
                          data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
bwqs.masld.pgs.1c_fit <- bwqs.masld.pgs.1c$summary_fit
## model output 
bwqs.masld.pgs.1c_fit <- bwqs.masld.pgs.1c$summary_fit
save(bwqs.masld.pgs.1c_fit, file = "pgs/bwqs.masld.pgs.1c.Rdata")
write.xlsx(bwqs.masld.pgs.1c_fit, "pgs/bwqs.masld.pgs.1c.xlsx")

# Adolescence Mix 1
bwqs.masld.pgs.1a.1 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                              bw + bf + fishadol_tert + ethn2c,
                            mix_name = adolmix1,
                            data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.pgs.1a.1_fit <- bwqs.masld.pgs.1a.1$summary_fit
save(bwqs.masld.pgs.1a.1_fit, file = "pgs/bwqs.masld.pgs.1a.1.Rdata")
write.xlsx(bwqs.masld.pgs.1a.1_fit, "pgs/bwqs.masld.pgs.1a.1.xlsx")

# Adolescence Mix 2
bwqs.masld.pgs.1a.2 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                              bw + bf + fishadol_tert + ethn2c,
                            mix_name = adolmix2,
                            data = data1.pgs, q=4, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.masld.pgs.1a.2_fit <- bwqs.masld.pgs.1a.2$summary_fit
save(bwqs.masld.pgs.1a.2_fit, file = "pgs/bwqs.masld.pgs.1a.2.Rdata")
write.xlsx(bwqs.masld.pgs.1a.2_fit, "pgs/bwqs.masld.pgs.1a.2.xlsx")



##############################################
# SENSITIVITY MODEL #4 (BWQS in tertiles instead of quantiles) ####
# 1. ALT ####

# Prenatal
bwqs.tert.alt.1m <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                           native + fishpreg_tert,
                         mix_name = premix,
                         data = data1, q=3, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output 
bwqs.tert.alt.1m_fit <- bwqs.tert.alt.1m$summary_fit
save(bwqs.tert.alt.1m_fit, file = "tertiles/bwqs.tert.alt.1m.Rdata")
write.xlsx(bwqs.tert.alt.1m_fit, "tertiles/bwqs.tert.alt.1m.xlsx")



# Childhood
bwqs.tert.alt.1c <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                           bw + bf + fishchild_tert + ethn2c,
                         mix_name = childmix,
                         data = data1, q=3, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
bwqs.tert.alt.1c_fit <- bwqs.tert.alt.1c$summary_fit
## model output 
bwqs.tert.alt.1c_fit <- bwqs.tert.alt.1c$summary_fit
save(bwqs.tert.alt.1c_fit, file = "tertiles/bwqs.tert.alt.1c.Rdata")
write.xlsx(bwqs.tert.alt.1c_fit, "tertiles/bwqs.tert.alt.1c.xlsx")



# Adolescence Mix 1
bwqs.tert.alt.1a.1 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             bw + bf + fishadol_tert + ethn2c,
                           mix_name = adolmix1,
                           data = data1, q=3, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.tert.alt.1a.1_fit <- bwqs.tert.alt.1a.1$summary_fit
save(bwqs.tert.alt.1a.1_fit, file = "tertiles/bwqs.tert.alt.1a.1.Rdata")
write.xlsx(bwqs.tert.alt.1a.1_fit, "tertiles/bwqs.tert.alt.1a.1.xlsx")



# Adolescence Mix 2
bwqs.tert.alt.1a.2 <- bwqs(alt ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             bw + bf + fishadol_tert + ethn2c,
                           mix_name = adolmix1,
                           data = data1, q=3, c_int = c(0.025, 0.975), family = "gaussian", seed = 1234)
## model output
bwqs.tert.alt.1a.2_fit <- bwqs.tert.alt.1a.2$summary_fit
save(bwqs.tert.alt.1a.2_fit, file = "tertiles/bwqs.tert.alt.1a.2.Rdata")
write.xlsx(bwqs.tert.alt.1a.2_fit, "tertiles/bwqs.tert.alt.1a.2.xlsx")

# 2. MASLD ####

# Prenatal
bwqs.tert.masld.1m <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             native + fishpreg_tert,
                           mix_name = premix,
                           data = data1, q=3, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output 
bwqs.tert.masld.1m_fit <- bwqs.tert.masld.1m$summary_fit
save(bwqs.tert.masld.1m_fit, file = "tertiles/bwqs.tert.masld.1m.Rdata")
write.xlsx(bwqs.tert.masld.1m_fit, "tertiles/bwqs.tert.masld.1m.xlsx")

# Childhood
bwqs.tert.masld.1c <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                             bw + bf + fishchild_tert + ethn2c,
                           mix_name = childmix,
                           data = data1, q=3, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
bwqs.tert.masld.1c_fit <- bwqs.tert.masld.1c$summary_fit
## model output 
bwqs.tert.masld.1c_fit <- bwqs.tert.masld.1c$summary_fit
save(bwqs.tert.masld.1c_fit, file = "tertiles/bwqs.tert.masld.1c.Rdata")
write.xlsx(bwqs.tert.masld.1c_fit, "tertiles/bwqs.tert.masld.1c.xlsx")

# Adolescence Mix 1
bwqs.tert.masld.1a.1 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                               bw + bf + fishadol_tert + ethn2c,
                             mix_name = adolmix1,
                             data = data1, q=3, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.tert.masld.1a.1_fit <- bwqs.tert.masld.1a.1$summary_fit
save(bwqs.tert.masld.1a.1_fit, file = "tertiles/bwqs.tert.masld.1a.1.Rdata")
write.xlsx(bwqs.tert.masld.1a.1_fit, "tertiles/bwqs.tert.masld.1a.1.xlsx")

# Adolescence Mix 2
bwqs.tert.masld.1a.2 <- bwqs(masld_2c ~ cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                               bw + bf + fishadol_tert + ethn2c,
                             mix_name = adolmix2,
                             data = data1, q=3, c_int = c(0.025, 0.975), family = "binomial", seed = 1234)
## model output
bwqs.tert.masld.1a.2_fit <- bwqs.tert.masld.1a.2$summary_fit
save(bwqs.tert.masld.1a.2_fit, file = "tertiles/bwqs.tert.masld.1a.2.Rdata")
write.xlsx(bwqs.tert.masld.1a.2_fit, "tertiles/bwqs.tert.masld.1a.2.xlsx")

