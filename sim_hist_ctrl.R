library(dplyr)
library(rstan)

S = 5
j = 1
Ns  = c(1000, 1000, 1000, 1000, 1000)
N  = sum(Ns)
theta = 2
delta = c(0, 0, 0, 1, -.5)
beta = .1
mu_x = c(10,10, 10, 15, 15)
sig = .1
y = NULL
X = NULL
Z = NULL
study = NULL
for (s in 1:S) {
  if (s == S | s==2) {
    x = c(rnorm(floor(Ns[s]/3), mu_x[4], 1), rnorm(Ns[s] - floor(Ns[s]/3), mu_x[1], 1))
  } else x = rnorm(Ns[s], mu_x[s], 1)
  z = rbinom(Ns[s], 1, .5)
  Z = c(Z, z)
  y = c(y, beta*x + theta*z + delta[s])
  X = c(X, x)
  study = c(study, rep(s, Ns[s]))
}
y = y + rnorm(length(y), 0, sig)
study = as.factor(study)
STUDY = model.matrix(~0+study)
df = data.frame(cbind(study, STUDY, X, y, Z))


dat = cbind(X, y)[Z==0,]
X1 = X[study == j & Z == 0]
y1 = y[study == j & Z == 0]
dat1 = cbind(X1, y1)
mu = apply(dat1, 2, mean)
cov = var(dat1)
d = apply(dat, 1, mahalanobis, center = mu, cov = cov)
w = 1-(d-min(d))/(max(d)-min(d))


X2 = X[study == j | Z == 0]
y2 = y[study == j | Z == 0]
Z2 = Z[study == j | Z == 0]
STUDY2 = STUDY[study == j | Z == 0,]
data = list(N = length(y2), Nw = length(w), y = y2, x = X2, z = Z2, w = w)
fit = stan(file='GaussCT_wWeights_wCovariates.stan', data=data, chains=0)
data = list(N = length(y2), S = S, y = y2, x = X2, z = Z2,  Z = STUDY2)
fit0 = stan(file='dbCT_Gauss_wCovariates.stan', data=data, chains=0)
###########################################################
df_sim = NULL
for (k in 1:100) {
  y = NULL
  X = NULL
  Z = NULL
  study = NULL
  for (s in 1:S) {
    if (s == S | s==2) {
      x = c(rnorm(floor(Ns[s]/3), mu_x[4], 1), rnorm(Ns[s] - floor(Ns[s]/3), mu_x[1], 1))
    } else x = rnorm(Ns[s], mu_x[s], 1)
    z = rbinom(Ns[s], 1, .5)
    Z = c(Z, z)
    y = c(y, beta*x + theta*z + delta[s])
    X = c(X, x)
    study = c(study, rep(s, Ns[s]))
  }
  y = y + rnorm(length(y), 0, sig)
  study = as.factor(study)
  STUDY = model.matrix(~0+study)
  df = data.frame(cbind(study, STUDY, X, y, Z))
  
  dat = cbind(X, y)
  X1 = X[study == j & Z == 0]
  y1 = y[study == j & Z == 0]
  dat1 = cbind(X1, y1)
  mu = apply(dat1, 2, mean)
  cov = var(dat1)
  d = apply(dat, 1, mahalanobis, center = mu, cov = cov)
  w = 1-(d-min(d))/(max(d)-min(d))
  w[Z==1] = 1
  w1 = w[study==j & Z==0]
  w[study==j] = 1
  eps = 1- quantile(w1, p = 0.05)
  
  X2 = X[study == j | Z == 0]
  y2 = y[study == j | Z == 0]
  Z2 = Z[study == j | Z == 0]
  w02 = w[study == j | Z == 0]
  study2 = study[study == j | Z == 0]
  STUDY2 = STUDY[study == j | Z == 0,]
  N2 = length(y2)
  data1 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w02)
  fit1 = stan('fit'=fit, 'data'=data1, warmup=1000, iter=2000, chains=2)
  fit_ss1 = extract(fit1)
  theta1 = fit_ss1$theta
  
  w0 = rep(1, N2)
  data2 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w0)
  fit2 = stan('fit'=fit, 'data'=data2, warmup=1000, iter=2000, chains=2)
  fit_ss2 = extract(fit2)
  theta2 = fit_ss2$theta
  
  w2 = w02
  w2[w2<(1-eps)] = 0
  data4 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w2)
  fit4 = stan('fit'=fit, 'data'=data4, warmup=1000, iter=2000, chains=2)
  fit_ss4 = extract(fit4)
  theta4 = fit_ss4$theta
  
  w3 = w02
  for (s in 1:S) w3[study2 == s] = mean(w3[study2 == s])
  data6 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w3)
  fit6 = stan('fit'=fit, 'data'=data6, warmup=1000, iter=2000, chains=2)
  fit_ss6 = extract(fit6)
  theta6 = fit_ss6$theta
  
  w4 = w3
  w4[w4<(1-eps)] = 0
  data7 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w4)
  fit7 = stan('fit'=fit, 'data'=data7, warmup=1000, iter=2000, chains=2)
  fit_ss7 = extract(fit7)
  theta7 = fit_ss7$theta
  
  y1 = y[study2 == j]
  w1 = w[study2 == j]
  X1 = X[study2 == j]
  Z1 = Z[study2 == j]
  N1 = length(y1)
  data3 = list(N = N1, Nw = length(w1), y = y1, x = X1, z = Z1, w = w1)
  fit3 = stan('fit'=fit, 'data'=data3, warmup=1000, iter=2000, chains=2)
  fit_ss3 = extract(fit3)
  theta3 = fit_ss3$theta
  
  
  data5 = list(N = N2, S = S, y = y2, x = X2, z = Z2,  Z = STUDY2)
  fit5 = stan('fit'=fit0, 'data'=data5, warmup=1000, iter=2000, chains=2)
  fit_ss5 = extract(fit5)
  theta5 = fit_ss5$theta
  df = data.frame(theta = c(theta1, theta2, theta3, theta4, theta5, theta6, theta7), 
                  method = c(rep('weighted',2000),rep( 'full', 2000), 
                             rep('none', 2000), rep('truncated w', 2000),
                             rep('db', 2000), rep('PowerPrior', 2000),
                             rep('truncated pp', 2000)))
  df_summary = df %>%
    group_by(method) %>%
    dplyr::summarize(meanTheta = mean(theta), lowTheta = quantile(theta, 0.025), 
                     uppTheta = quantile(theta, 0.975), RMSE = sqrt(mean((theta - 2)^2)))
  sim = rep(k, nrow(df_summary))
  df_sim = data.frame(rbind(df_sim, cbind(sim, df_summary)))
}

save(df_sim, file = 'sim_mahala_HC.Rdata')
load('sim_mahala_HC.Rdata')
df_sim$method = plyr::mapvalues(df_sim$method, c("db","full","none","PowerPrior",
                                                 "truncated pp", "truncated w","weighted" ), 
                                c('CP', 'FH', 'NP', 'PP', 'TPP', 'TIW', 'IW'))
df_sim$method = factor(df_sim$method, levels(df_sim$method)[c(2, 7, 4, 1, 3, 5, 6)])
ggplot(df_sim, aes(x = method, y = RMSE, fill = method)) + 
  geom_boxplot(size = .35, width = .5, outlier.alpha = 0, color = 'darkgrey') +
  scale_fill_brewer(palette = 'Set2') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
ggsave(file = 'sim_all_2.png')

df <- data.frame(x=1:8, y=1, col=letters[1:8])
g <- ggplot(df, aes(x=x, y=y, color=col)) + geom_point(size=5) +
  scale_color_brewer(palette="Set2")
set2 <- ggplot_build(g)$data[[1]]$colour

df_sim0 = df_sim %>%
  filter(method %in% c('NP', 'CP', 'TPP', 'TIW'))
ggplot(df_sim0, aes(x = method, y = RMSE, fill = method)) + 
  geom_boxplot(size = .5, width = .25, outlier.alpha = 0, color = 'darkgrey') +
  #ylim(c(.025,.1)) +
  scale_fill_manual(values = set2[4:7]) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
ggsave(file = 'sim_some_2.png')

################bias and variance####################
df1 = df_sim %>%
  mutate(bias = (meanTheta - 2), IQE = (uppTheta - lowTheta))
df_sim$method = plyr::mapvalues(df_sim$method, c("db","full","none","PowerPrior",
                                                 "truncated pp", "truncated w","weighted" ), 
                                c('CP', 'FH', 'NP', 'PP', 'TPP', 'TIW', 'IW'))
df1$method = factor(df1$method, levels(df_sim$method)[c(2, 7, 4, 1, 3, 5, 6)])
ggplot(df1, aes(x = method, y = bias, fill = method)) + 
  geom_boxplot(size = .35, width = .5, outlier.alpha = 0, color = 'darkgrey') +
  scale_fill_brewer(palette = 'Set2') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
ggsave(file = 'sim_all_bias.png')

ggplot(df1, aes(x = method, y = IQE, fill = method)) + 
  geom_boxplot(size = .35, width = .5, outlier.alpha = 0, color = 'darkgrey') +
  scale_fill_brewer(palette = 'Set2') +
  ylab('95% credible interval length') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
ggsave(file = 'sim_all_IQE.png')

