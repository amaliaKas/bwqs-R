##############################################
##########Project: HELIX study - adolescence - PFAS mixtures and liver outcomes
##########Creator: Amalia Kasmiridou
##########Aim: BWQS effect modification Analyses for imp=1
##########Last update: 22/07/2026 
##############################################
# Set working directory & load the dataset ####
setwd("~/BWQS/")
load("data_pfas_liver.Rdata")

# Imputed dataset #1
data1 <- subset(data_pfas_liver, .imp==1)

##############################################
# EFFECT MODIFICATION ####
##############################################
# Prepare data ####
d <- as.data.frame(data1)

d$sex <- as.numeric(ifelse(d$sex == "male", 1, 0)) #females are 0, males are 1
d$upf_2c <- as.numeric(ifelse(d$upf_2c == 1, 0, 1)) #low/mid are 0, high are 1
d$pgs0027_2c <- as.numeric(ifelse(d$pgs0027_2c == 1, 0, 1)) #low/mid are 0, high are 1
d$pgs2071_2c <- as.numeric(ifelse(d$pgs2071_2c == 1, 0, 1)) #low/mid are 0, high are 1
d$alt <- as.numeric(d$alt)
d$hs2_masld_2c <- as.numeric(d$masld_2c) #low risk are 0, high risk ("cases") are 1

##############################################
# MODEL ####


library('rstan')
library('ggplot2')
library('BWQS')
library('mvtnorm')
library('tidyverse')
library('DT')

model_bwqs_gaussian_group_lasso <- "data {

int<lower=0> N;              // Sample size
int<lower=0> C1;             // number of element of first mix
int<lower=0> K;              // number of covariates
matrix[N,C1] XC1;	           // matrix of first mix
vector[C1] DalpC1;           // vector of the Dirichlet coefficients for first mix
matrix[N,K] KV;	             // matrix of covariates
vector[N] interaction_term;         // interaction variable
real y[N];                   // outcome gaussian variable
}

parameters {

real <lower=0> sigma;
real mu;               // intercepts
real beta;             // overall mixture coeffs
real beta_int;         // main interaction coeff with mixture index
real delta_int;        // coeffs for interaction term 
vector[K] delta;       // covariates coefficients
simplex[C1] WC1;       // weights of first mix
}

transformed parameters {

vector[N] Xb;
real bbeta0;  // baseline zero group
real bbeta1;  // group above
Xb = mu + (XC1*WC1)*beta + interaction_term*delta_int  +  ((XC1*WC1).*interaction_term)*beta_int + KV*delta;
bbeta0 = beta;
bbeta1 = beta + beta_int;
}

model {

mu ~ normal(0, 100);
sigma ~ inv_gamma(0.01,0.01);
beta ~ normal(0,100);
beta_int ~ normal(0,100);
delta_int ~ normal(0,100);
delta ~ normal(0,100);

WC1 ~ dirichlet(DalpC1);
y ~  normal(Xb, sigma);
}
"
m_lasso <- rstan::stan_model(model_code =  model_bwqs_gaussian_group_lasso)

##############################################
# 1. SEX ####
# PRENATAL ####

data <- d[,c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'native','fishpreg_tert',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       native + fishpreg_tert)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "sex"

mix_name_1 <- c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m')

X1 = BWQS::quantile_split(data=data, mix_name=mix_name_1, q)[,mix_name_1]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_1),
  XC1 = cbind(X1),
  DalpC1 = rep(1, length(mix_name_1)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.sex.1m.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.sex.1m.xlsx")

# CHILDHOOD ####

data <- d[,c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishchild_tert', 'ethn2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishchild_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "sex"

mix_name_2 <- c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c')

X2 = BWQS::quantile_split(data=data, mix_name=mix_name_2, q)[,mix_name_2]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_2),
  XC1 = cbind(X2),
  DalpC1 = rep(1, length(mix_name_2)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]', 'WC1[5]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.sex.1c.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.sex.1c.xlsx")

# ADOLESCENCE 1 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "sex"

mix_name_3 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a')

X3 = BWQS::quantile_split(data=data, mix_name=mix_name_3, q)[,mix_name_3]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_3),
  XC1 = cbind(X3),
  DalpC1 = rep(1, length(mix_name_3)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

# rownames(sum_fit_lasso) <- c("Intercept","Chemical Mixture")
datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.sex.1a.1.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.sex.1a.1.xlsx")

# ADOLESCENCE 2 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "sex"

mix_name_4 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
                'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a')

X4 = BWQS::quantile_split(data=data, mix_name=mix_name_4, q)[,mix_name_4]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_4),
  XC1 = cbind(X4),
  DalpC1 = rep(1, length(mix_name_4)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]',
                'WC1[5]', 'WC1[6]', 'WC1[7]', 'WC1[8]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

# rownames(sum_fit_lasso) <- c("Intercept","Chemical Mixture")
datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.sex.1a.2.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.sex.1a.2.xlsx")

##############################################
# 2. UPF ####
# PRENATAL ####

data <- d[,c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'native','fishpreg_tert', 
             # plus the interaction variables
             'upf_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       native + fishpreg_tert)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "upf_2c"

mix_name_1 <- c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m')

X1 = BWQS::quantile_split(data=data, mix_name=mix_name_1, q)[,mix_name_1]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_1),
  XC1 = cbind(X1),
  DalpC1 = rep(1, length(mix_name_1)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.upf.1m.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.upf.1m.xlsx")


# CHILDHOOD ####

data <- d[,c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishchild_tert', 'ethn2c',
             # plus the interaction variables
             'upf_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishchild_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "upf_2c"

mix_name_2 <- c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c')

X2 = BWQS::quantile_split(data=data, mix_name=mix_name_2, q)[,mix_name_2]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_2),
  XC1 = cbind(X2),
  DalpC1 = rep(1, length(mix_name_2)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)


fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]', 'WC1[5]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.upf.1c.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.upf.1c.xlsx")

# ADOLESCENCE 1 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             # plus the interaction variables
             'upf_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "upf_2c"

mix_name_3 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a')

X3 = BWQS::quantile_split(data=data, mix_name=mix_name_3, q)[,mix_name_3]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_3),
  XC1 = cbind(X3),
  DalpC1 = rep(1, length(mix_name_3)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.upf.1a.1.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.upf.1a.1.xlsx")

# ADOLESCENCE 2 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             # plus the interaction variables
             'upf_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "upf_2c"

mix_name_4 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
                'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a')

X4 = BWQS::quantile_split(data=data, mix_name=mix_name_4, q)[,mix_name_4]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_4),
  XC1 = cbind(X4),
  DalpC1 = rep(1, length(mix_name_4)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]',
                'WC1[5]', 'WC1[6]', 'WC1[7]', 'WC1[8]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.upf.1a.2.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.upf.1a.2.xlsx")


##############################################
# 3. PGS BMI ####
# PRENATAL ####

data <- d[,c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'native','fishpreg_tert', 
             # plus the interaction variables
             'pgs0027_2c', 
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       native + fishpreg_tert)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs0027_2c"

mix_name_1 <- c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m')

X1 = BWQS::quantile_split(data=data, mix_name=mix_name_1, q)[,mix_name_1]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_1),
  XC1 = cbind(X1),
  DalpC1 = rep(1, length(mix_name_1)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.bmi.1m.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.bmi.1m.xlsx")


# CHILDHOOD ####

data <- d[,c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishchild_tert', 'ethn2c',
             # plus the interaction variables
             'pgs0027_2c', 
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishchild_tert + ethn2c)


y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs0027_2c"

mix_name_2 <- c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c')

X2 = BWQS::quantile_split(data=data, mix_name=mix_name_2, q)[,mix_name_2]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_2),
  XC1 = cbind(X2),
  DalpC1 = rep(1, length(mix_name_2)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]', 'WC1[5]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.bmi.1c.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.bmi.1c.xlsx")

# ADOLESCENCE 1 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             # plus the interaction variables
             'pgs0027_2c', 
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs0027_2c"

mix_name_3 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a')

X3 = BWQS::quantile_split(data=data, mix_name=mix_name_3, q)[,mix_name_3]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_3),
  XC1 = cbind(X3),
  DalpC1 = rep(1, length(mix_name_3)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.bmi.1a.1.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.bmi.1a.1.xlsx")

# ADOLESCENCE 2 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             # plus the interaction variables
             'pgs0027_2c', 
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs0027_2c"

mix_name_4 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
                'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a')

X4 = BWQS::quantile_split(data=data, mix_name=mix_name_4, q)[,mix_name_4]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_4),
  XC1 = cbind(X4),
  DalpC1 = rep(1, length(mix_name_4)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 2,
                             iter = 1e3,
                             thin = 1,
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]',
                'WC1[5]', 'WC1[6]', 'WC1[7]', 'WC1[8]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.bmi.1a.2.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.bmi.1a.2.xlsx")


##############################################
# 4. PGS LD (CHANGED PARAMETERS) ####
# PRENATAL ####

data <- d[,c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'native','fishpreg_tert', 
             # plus the interaction variables
             'pgs2071_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       native + fishpreg_tert)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs2071_2c"

mix_name_1 <- c('log2_pfhxs_m','log2_pfna_m','log2_pfoa_m','log2_pfos_m')

X1 = BWQS::quantile_split(data=data, mix_name=mix_name_1, q)[,mix_name_1]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_1),
  XC1 = cbind(X1),
  DalpC1 = rep(1, length(mix_name_1)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 3,        # number of simultaneous processes
                             iter = 1e4,        # n of times 
                             thin = 5,          # we chose each 5th iter estimates 
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.ld.2.1m.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.ld.2.1m.xlsx")


# CHILDHOOD ####

data <- d[,c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishchild_tert', 'ethn2c',
             # plus the interaction variables
             'pgs2071_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishchild_tert + ethn2c)


y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs2071_2c"

mix_name_2 <- c('log2_pfhxs_c','log2_pfna_c','log2_pfoa_c','log2_pfos_c', 'log2_pfunda_c')

X2 = BWQS::quantile_split(data=data, mix_name=mix_name_2, q)[,mix_name_2]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_2),
  XC1 = cbind(X2),
  DalpC1 = rep(1, length(mix_name_2)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 3,        # number of simultaneous processes
                             iter = 1e4,        # n of times 
                             thin = 5,          # we chose each 5th iter estimates 
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]', 'WC1[5]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.ld.2.1c.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.ld.2.1c.xlsx")

# ADOLESCENCE 1 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             # plus the interaction variables
             'pgs2071_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)

y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs2071_2c"

mix_name_3 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a')

X3 = BWQS::quantile_split(data=data, mix_name=mix_name_3, q)[,mix_name_3]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_3),
  XC1 = cbind(X3),
  DalpC1 = rep(1, length(mix_name_3)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 3,        # number of simultaneous processes
                             iter = 1e4,        # n of times 
                             thin = 5,          # we chose each 5th iter estimates 
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.ld.2.1a.1.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.ld.2.1a.1.xlsx")

# ADOLESCENCE 2 ####

data <- d[,c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
             'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a',
             'cohort', 'sex', 'mbmi','mage','medu','parity','msmok','adolage', 'ipw',
             'bw', 'bf', 'fishadol_tert', 'ethn2c',
             # plus the interaction variables
             'pgs2071_2c',
             "alt")]

data <- na.omit(data)

q = 4  # number of quantiles

# Change covariates here
formula = as.formula(alt ~ # write only the name of the covariates
                       cohort + sex + mbmi + mage + medu + parity + msmok + adolage + ipw +
                       bw + bf + fishadol_tert + ethn2c)


y_name  <- all.vars(formula)[1]
KV_name <- all.vars(formula)[-1]

interaction.name <- "pgs2071_2c"

mix_name_4 <- c('log2_pfhxs_a','log2_pfna_a','log2_pfoa_a','log2_pfos_a',
                'log2_pfunda_a', 'log2_pfhpa_a', 'log2_pfda_a', 'log2_x9clpfesa_a')

X4 = BWQS::quantile_split(data=data, mix_name=mix_name_4, q)[,mix_name_4]

# run the fit_lasso model

data_reg <- list(
  
  N   = nrow(data),
  
  C1  = length(mix_name_4),
  XC1 = cbind(X4),
  DalpC1 = rep(1, length(mix_name_4)),
  
  interaction_term = data[,interaction.name],
  
  KV = data[,KV_name],
  K   = length(KV_name),
  
  y = as.vector(data[,y_name])
)

fit_lasso <- rstan::sampling(m_lasso,
                             data = data_reg,
                             chains = 3,        # number of simultaneous processes
                             iter = 1e4,        # n of times 
                             thin = 5,          # we chose each 5th iter estimates 
                             refresh = 0, verbose = T,
                             control=list(max_treedepth = 20,
                                          adapt_delta = 0.999999999999999))

s <- summary(fit_lasso)
s$summary

parameters <- c('mu','beta', 'bbeta1', 'beta_int', 'delta_int', 'WC1[1]', 'WC1[2]', 'WC1[3]', 'WC1[4]',
                'WC1[5]', 'WC1[6]', 'WC1[7]', 'WC1[8]')
sum_fit_lasso <- round(summary(fit_lasso,pars = parameters,
                               probs = c(0.025, 0.975))$summary,5)

datatable(sum_fit_lasso, class = 'cell-border stripe',options = list(dom = 't', pageLength = 20))

save(sum_fit_lasso, file = "eff-mod/bwqs.alt.ld.2.1a.2.Rdata")
write.xlsx(sum_fit_lasso, "eff-mod/bwqs.alt.ld.2.1a.2.xlsx")

