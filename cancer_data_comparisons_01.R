library(rstan)
library(dplyr)
library(ggplot2)


load('TTE_4study.Rdata')
j = 4
Z = tte_dat[,c('STUDY1', 'STUDY2','STUDY3','STUDY4')]
S = ncol(Z)
X = tte_dat[,c('AGE', 'SEX', 'RACE1', 'RACE2','RACE3', 'RACE4')]
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
X = tte_ctrl[,c('AGE', 'SEX', 'RACE1', 'RACE2','RACE3', 'RACE4')]
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
mu = apply(dat1, 2,mean)
cov = var(dat1)
dat = cbind(X, y)

####### ML weights ##########
x1 = X1[,1]
x2 = X1[,2:5]


data001 = list(C = 5, N = length(y1), y = y1, nu = nu1, x1 = x1, x2 = x2)
fit001 = stan(file='Gauss_similarity.stan', data=data001, chains=0)
fit001 = stan('fit'=fit001, 'data'=data001, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit001)


weight_fun = function(y, nu, x1, x2, pars) {
  theta1 = pars$theta1
  theta2 = pars$theta2
  beta1 = pars$beta1
  beta2 = pars$beta2
  mu_x = pars$mu_x
  p_x = pars$p_x
  sig_x = pars$sig_x
  sig_y = pars$sig_y
  p_nu = pars$p_nu
  lambda = beta1*x1 + c(beta2%*%t(x2))
  llx = matrix(0, length(theta1), length(x2))
  for (n in 1:length(theta)) {
    for (k in 1:length(x2)) llx[n,k] = dbinom(as.numeric(x2[k]), 1, p_x[n,k], log = T)
  }
  if (nu == 1) {
    lly = dnorm(y, lambda+theta1, sig_y, log = T)
  } else lly = dnorm(y, lambda+theta2, sig_y, log = T)
  ll = lly + dnorm(x1, mu_x, sig_x, log = T) + apply(llx, 1, sum) + dbinom(nu, 1, p_nu, log = T)
  return(sum(exp(ll)))
}

w0 = c()
for (n in 1:length(y)) w0[n] = weight_fun(y[n], nu[n], X[n,1], X[n,2:5], pars = fit_ss)
#w0[nu==0] = 0
w2 = (w0-min(w0))/(max(w0)-min(w0))
eps2 = quantile(w2[study==j], p = .05)
w2[study==j] = 1


data = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w2, y_t = y_t, nu = nu, nu_t = nu_t, 
            X = X, X_t = X_t, Z = Z)
fit1 = stan(file='survival_wCovariates_wWeights.stan', data=data, chains=0)
data = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, y_t = y_t, nu = nu, nu_t = nu_t, 
            X = X, X_t = X_t, Z = Z)
fit2 = stan(file='db_survival_wCovariates_01.stan', data=data, chains=0)

############################################################################################
### Comensurate Prior
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

#############################################################################################
### No Prior

w_NP = w2
w_NP[df$STUDY4 == 0] = 0
data3 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_NP, y_t = y_t, nu = nu, nu_t = nu_t, 
             X = X, X_t = X_t, Z = Z)
fit = stan('fit'=fit1, 'data'=data3, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit)
HR_NP = exp(fit_ss$delta)
theta_NP = fit_ss$theta
alpha_NP = fit_ss$alpha
beta_NP = fit_ss$beta
#############################################################################################
### Power Prior

w_PP = w2
w_PP[df$STUDY1 == 1] = mean(w2[df$STUDY1 == 1])
w_PP[df$STUDY2 == 1] = mean(w2[df$STUDY2 == 1])
w_PP[df$STUDY3 == 1] = mean(w2[df$STUDY3 == 1])
w_PP[df$STUDY4 == 1] = mean(w2[df$STUDY4 == 1])
data4 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_PP, y_t = y_t, nu = nu, nu_t = nu_t,
             X = X, X_t = X_t, Z = Z)
fit = stan('fit'=fit1, 'data'=data4, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit)
HR_PP = exp(fit_ss$delta)

#############################################################################################
### Truncated Power Prior

# w_TPP = w_PP
# w_TPP[w_TPP<thresh] = 0
# data5 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_TPP, y_t = y_t, nu = nu, nu_t = nu_t, 
#              X = X, X_t = X_t, Z = Z)
# fit = stan('fit'=fit1, 'data'=data5, warmup=1000, iter=2000, chains=2)
# fit_ss = extract(fit)
# HR_TPP = exp(fit_ss$delta)
# theta_TPP = fit_ss$theta
# alpha_TPP = fit_ss$alpha
# beta_TPP = fit_ss$beta
#############################################################################################
### Individually Weighted Prior

 w_IW = w2
# data6 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_IW, y_t = y_t, nu = nu, nu_t = nu_t, 
#              X = X, X_t = X_t, Z = Z)
# fit = stan('fit'=fit1, 'data'=data6, warmup=1000, iter=2000, chains=2)
# fit_ss = extract(fit)
# HR_IW = exp(fit_ss$delta)

#############################################################################################
### Truncated Individually Weighted Prior

w_TIW = w_IW
w_TIW[w_TIW<thresh] = 0
data7 = list(j = j, C = C, S = S, N = N, Nt = Nt, y = y, w = w_TIW, y_t = y_t, nu = nu, nu_t = nu_t, 
             X = X, X_t = X_t, Z = Z)
fit = stan('fit'=fit1, 'data'=data7, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit)
HR_TIW = exp(fit_ss$delta)
theta_TIW = fit_ss$theta
alpha_TIW = fit_ss$alpha
beta_TIW = fit_ss$beta
############################################################################################
# df1 = data.frame(HR = c(HR_CP, HR_FH, HR_NP, HR_PP, HR_TPP, HR_IW, HR_TIW), 
#                 method = c(rep('CP',2000), rep( 'FH', 2000), rep('NP', 2000), 
#                            rep('PP', 2000), rep('TPP', 2000), rep('IW', 2000),
#                            rep('TIW', 2000)))

df = data.frame(HR = c(HR_NP, HR_TIW), 
                   theta = c(theta_NP, theta_TIW),
                alpha = c(alpha_NP, alpha_TIW),
                 method = c(rep('NP', 2000), 
                            rep('TIW', 2000)))

df2 = data.frame(HR = c(HR_CP, HR_FH, HR_NP, HR_PP,  HR_TIW),
                method = c(rep('CP',2000), rep( 'FH', 2000), rep('NP', 2000),
                           rep('PP', 2000), rep('TIW', 2000)))

df_summary2 = df2 %>%
  group_by(method) %>%
  dplyr::summarize(HRmean = mean(HR), HRlow = quantile(HR, 0.025), HRupp = quantile(HR, 0.975))


df_summary = df %>%
  group_by(method) %>%
  dplyr::summarize(HRmean = mean(HR), HRlow = quantile(HR, 0.025), HRupp = quantile(HR, 0.975),
                   theta_mean = mean(exp(theta)), theta_low = quantile(exp(theta), 0.025), 
                   theta_upp = quantile(exp(theta), 0.975), alpha_mean = mean(alpha), 
                   alpha_low = quantile(alpha, 0.025), 
                   alpha_upp = quantile(alpha, 0.975))

df_summary2$method = factor(df_summary2$method, levels(df_summary2$method)[c(3,5,1,4,2)])
ggplot(df_summary2, aes(x = HRmean, y = method, color = method)) + geom_point(size = 2) +
  geom_errorbarh(aes(xmax = HRupp, xmin = HRlow), height = 0, size = 1) + 
  xlab('HR') + 
  scale_color_brewer(palette = 'Set1') +
  theme(axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = 'bold'),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = 'bold'),
    legend.position = 'none')
ggsave(file = 'cancer_data_HR2.png')

ggplot(df_summary, aes(x = log(theta_mean), y = method, color = method)) + geom_point(size = 2) +
  geom_errorbarh(aes(xmax = log(theta_upp), xmin = log(theta_low)), height = 0, size = 1) + 
  xlab(expression(theta)) + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = 'bold'),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = 'bold'),
    legend.position = 'none')
ggsave(file = 'cancer_data_theta.png',width = 7, height = 3)

ggplot(df_summary, aes(x = alpha_mean, y = method, color = method)) + geom_point(size = 2) +
  geom_errorbarh(aes(xmax = alpha_upp, xmin = alpha_low), height = 0, size = 1) + 
  xlab('shape parameter') + 
  scale_color_brewer(palette = 'Set1') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'none')
ggsave(file = 'cancer_data_alpha.png',width = 7, height = 3)

