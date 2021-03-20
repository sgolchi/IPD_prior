load('TTE_4study.Rdata')
library(dplyr)
library(rstan)
library(ggplot2)
library(mvtnorm)
j = 1
tte_dat$RACE14 = tte_dat$RACE1 | tte_dat$RACE4
Z = tte_dat[,c('STUDY1', 'STUDY2','STUDY3','STUDY4')]
S = ncol(Z)
X = tte_dat[,c('AGE', 'SEX', 'RACE2','RACE3', 'RACE14')]
A = ifelse(tte_dat[, 'arm'] == 'ctrl', 0, 1)
C = ncol(X)
y = tte_dat$SURV
N = length(y)

tte_ctrl = tte_dat %>%
  filter(arm == 'ctrl')
tte_trt = tte_dat %>%
  filter(arm == 'trt' & STUDY == j)
Z = tte_ctrl[,c('STUDY1', 'STUDY2','STUDY3','STUDY4')]
S = ncol(Z)
X = tte_ctrl[,c('AGE', 'SEX', 'RACE2','RACE3', 'RACE14')]
C = ncol(X)
y = tte_ctrl$SURV
nu = tte_ctrl$CENS
N = length(y)
y_t = tte_trt$SURV
nu_t = tte_trt$CENS
Nt = length(y_t)
X_c = X[which(tte_ctrl$STUDY == j),]
X_t = apply(X_c, 2, sample, size = Nt, replace = T)


STUDY = Z
df = data.frame(cbind(X, STUDY, y))

study = as.factor(tte_ctrl$STUDY)
X1 = X[study == j,]
y1 = y[study==j]
nu1 = nu[study==j]
dat1 = cbind(X1, y1)
mu = apply(dat1, 2, mean)
cov = var(dat1)
dat = cbind(X, y)

cat.v = cbind(tte_ctrl$SEX,  (tte_ctrl$RACE2 | tte_ctrl$RACE3 | tte_ctrl$RACE14), tte_ctrl$CENS)
cont.v = cbind(tte_ctrl$AGE, tte_ctrl$SURV)
out = fit_simmodel(cat.v, cont.v, which(study==j))
w = out$weights
#w01 = GMDw(cat.v, cont.v, s.ind = study==j) 
eps2 = quantile(w[study==j], p = .05)
df = data.frame(w = w, study = study)
df$study = plyr::mapvalues(df$study, c(1,2,3,4), c('STUDY57', 'INTEREST', 'PROCLAIM', 'ZODIAC'))
ggplot(df, aes(x = w)) + geom_histogram() + facet_grid(study ~.) + 
  geom_vline(xintercept = eps2, col = 'darkred')
ggsave(file = 'weights.png')

data = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w, y_t = y_t, nu = nu, nu_t = nu_t, 
            X = X, X_t = X_t, Z = Z)
fit1 = stan(file='survival_wCovariates_wWeights.stan', data=data, chains=0)
data = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, y_t = y_t, nu = nu, nu_t = nu_t, 
            X = X, X_t = X_t, Z = Z)
fit2 = stan(file='db_survival_wCovariates_01.stan', data=data, chains=0)
data = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, y_t = y_t, nu = nu, nu_t = nu_t, 
            X = X, X_t = X_t, Z = Z, s = as.integer(study))
fit3 = stan(file='survival_wCovariates_pp.stan', data=data, chains=0)

############################################################################################
### MAP
data1 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, y_t = y_t, nu = nu, nu_t = nu_t,
             X = X, X_t = X_t, Z = Z)
fit = stan('fit'=fit2, 'data'=data1, warmup=1000, iter=2000, chains=2)

fit_ss = extract(fit)
HR_CP = exp(fit_ss$delta)
theta_CP = fit_ss$theta[,4]
alpha_CP = fit_ss$alpha
beta_CP = fit_ss$beta
############################################################################################
### Full History
w_FH = rep(1,N)
data2 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_FH, y_t = y_t, nu = nu, nu_t = nu_t,
             X = X, X_t = X_t, Z = Z)
fit = stan('fit'=fit1, 'data'=data2, warmup=1000, iter=2000, chains=2)

fit_ss = extract(fit)
HR_FH = exp(fit_ss$delta)
delta_FH = fit_ss$delta
theta_FH = fit_ss$theta
alpha_FH = fit_ss$alpha
beta_FH = fit_ss$beta

#############################################################################################
### No Prior

w_NP = w
w_NP[study!=j] = 0
w_NP[study==j] = 1
data3 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_NP, y_t = y_t, nu = nu, nu_t = nu_t, 
             X = X, X_t = X_t, Z = Z)
fit_np = stan('fit'=fit1, 'data'=data3, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit_np)
HR_NP = exp(fit_ss$delta)
delta_Np = fit_ss$delta
theta_NP = fit_ss$theta
alpha_NP = fit_ss$alpha
beta_NP = fit_ss$beta

#############################################################################################
### Truncated Individually Weighted Prior

w_TIW = w
w_TIW[w_TIW<eps2] = 0
w_TIW[study==j] = 1
data7 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_TIW, y_t = y_t, nu = nu, nu_t = nu_t, 
             X = X, X_t = X_t, Z = Z)
fit_tiw = stan('fit'=fit1, 'data'=data7, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit_tiw)
HR_TIW = exp(fit_ss$delta)
delta_TIW = fit_ss$delta
theta_TIW = fit_ss$theta
alpha_TIW = fit_ss$alpha
beta_TIW = fit_ss$beta

#############################################################################################
### Power prior

data4 = list(j = 1, C = C, S = S, N = N, Nt = Nt, y = y, y_t = y_t, nu = nu, nu_t = nu_t, 
             X = X, X_t = X_t, Z = Z, s = as.integer(study))
fit_PP = stan('fit'=fit3, 'data'=data4, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit_PP)
HR_pp = exp(fit_ss$delta)
theta_pp = fit_ss$theta
alpha_pp = fit_ss$alpha
beta_pp = fit_ss$beta





df2 = data.frame(HR = c(HR_CP, HR_FH, HR_NP, HR_pp,  HR_TIW),
                 method = c(rep('CP',2000), rep( 'FH', 2000), rep('NP', 2000),
                            rep('PP', 2000), rep('TIW', 2000)))

df3 = data.frame(HR = c(HR_FH, HR_NP, HR_TIW, HR_CP),
                theta= c(theta_FH, theta_NP, theta_TIW, theta_CP),
                alpha = c(alpha_FH, alpha_NP, alpha_TIW, alpha_CP),
                 method = c(rep( 'FH', 2000), rep('NP', 2000),rep('TIW', 2000), 
                            rep('MAP', 2000)))


df_summary2 = df3 %>%
  group_by(method) %>%
  dplyr::summarize(HRmean = mean(HR), HRlow = quantile(HR, 0.025), HRupp = quantile(HR, 0.975),
                   thetamean = mean(theta), thetalow = quantile(theta, 0.025), thetaupp = quantile(theta, 0.975),
                   alphamean = mean(alpha), alphalow = quantile(alpha, 0.025), alphaupp = quantile(alpha, 0.975))

df_summary2$method = factor(df_summary2$method, levels = c('TIW','MAP', 'NP', 'FH'))


#plots
# ggplot(df_summary2, aes(x = HRmean, y = method)) + geom_point(size = 2, color = 'grey40') +
#   geom_errorbarh(aes(xmax = HRupp, xmin = HRlow), height = 0, size = 1, color = 'grey40') + 
#   xlab('HR') + 
#   scale_color_brewer(palette = 'Set1') +
#   theme_bw() +
#   geom_vline(xintercept = 1, color = 'darkred') + 
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'none')
# ggsave(file = 'CD_HRj1.png')
# 
# ggplot(df_summary2, aes(x = thetamean, y = method, color = method)) + geom_point(size = 2) +
#   geom_errorbarh(aes(xmax = thetaupp, xmin = thetalow), height = 0, size = 1) + 
#   xlab('theta') + 
#   scale_color_brewer(palette = 'Set1') +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'none')
# ggsave(file = 'CD_thetaj1.png')
# 
# ggplot(df_summary2, aes(x = alphamean, y = method, color = method)) + geom_point(size = 2) +
#   geom_errorbarh(aes(xmax = alphaupp, xmin = alphalow), height = 0, size = 1) + 
#   xlab('alpha') + 
#   scale_color_brewer(palette = 'Set1') +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'none')
# ggsave(file = 'CD_alphaj1.png')
# 
# 
# dfw = data.frame(weight = w, study = study)
# ggplot(dfw, aes(x = weight, fill = study)) + geom_histogram() + facet_grid(study~.) +
#   geom_vline(xintercept = eps2)

