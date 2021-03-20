library(dplyr)
library(rstan)

S = 2
j = 1
Ns  = c(100, 100)
N  = sum(Ns)
theta = 2
delta = c(0, 2)
beta = .2
mu_x = c(10, 10)
sig = 1
y = NULL
X = NULL
Z = NULL
study = NULL
for (s in 1:S) {
  if (s==2) {
    x = c(rnorm(floor(Ns[s]/2), mu_x[1], 1), rnorm(Ns[s] - floor(Ns[s]/2), mu_x[2], 1))
    del = c(rep(1, floor(Ns[s]/2)) * delta[1], rep(1, Ns[s] - floor(Ns[s]/2)) * delta[2])
  } else {
    x = rnorm(Ns[s], mu_x[s], 1)
    del = rep(1, Ns[s]) * delta[1]
  }
  z = rbinom(Ns[s], 1, .5)
  Z = c(Z, z)
  y = c(y, beta*x + theta*z + del)
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


X2 = X[study == j | Z == 0]
y2 = y[study == j | Z == 0]
Z2 = Z[study == j | Z == 0]
STUDY2 = STUDY[study == j | Z == 0,]
study2 = study[study == j | Z == 0]
w = w[study == j | Z == 0]
data = list(N = length(y2), Nw = length(w), y = y2, x = X2, z = Z2, w = w)
fit = stan(file='GaussCT_wWeights_wCovariates.stan', data=data, chains=0)
data = list(N = length(y2), S = S, y = y2, x = X2, z = Z2,  Z = STUDY2)
fit0 = stan(file='dbCT_Gauss_wCovariates.stan', data=data, chains=0)
data = list(N = length(y2), S = 5, s = as.integer(study2), y = y2, X = X2, z = Z2)
fitpp = stan(file='GaussCT_pp_wCovariates_1x.stan', data=data, chains=0)


bigG = function(a0, X0, y0, Xc, yc) {
  n = length(y)
  n0 = length(y0)
  p = ncol(X0)
  dmat = t(Xc)%*%Xc + a0*t(X0)%*%X0
  c = (n+a0*n0-p)*log(2*pi*sig^2) + log(det(dmat))
  b = t(yc)%*%yc + a0*t(y0)%*%y0
  a = t(Xc)%*%yc + a0*t(X0)%*%y0
  d = b - t(a)%*%solve(dmat)%*%a
  e = d/sig^2 + log(n0)/a0
  return(e)
}
###########################################################
sim_fun = function(delta_h1=3, delta_h2=3, mu_h1=1, mu_h2=1, delta_c = 0, mu_c = 1, 
                   theta_t = 0.5, beta = 1, sig = 1, Nc = 100, 
                   Nh = 100, M = 500) {
  df_sim = NULL
  for (k in 1:M) {
    y = NULL
    X = NULL
    Z = NULL
    study = NULL
    x = rnorm(Nc, mu_c, 1)
    del = rep(1, Nc) * delta_c
    z = rbinom(Nc, 1, .5)
    Z = c(Z, z)
    y = c(y, beta*x + theta_t*z + del)
    X = c(X, x)
    study = c(study, rep(1, Nc))
    x = c(rnorm(Nh - floor(Nh/4), mu_h1, 1), rnorm(floor(Nh/4), mu_h2, 1))
    del = c(rep(1, Nh - floor(Nh/4)) * delta_h1, rep(1, floor(Nh/4)) * delta_h2)
    z = rbinom(Nh, 1, .5)
    Z = c(Z, z)
    y = c(y, beta*x + theta_t*z + del)
    X = c(X, x)
    study = c(study, rep(2, Nh))
    y = y + rnorm(length(y), 0, sig)
    study = as.factor(study)
    STUDY = model.matrix(~0+study)
    df = data.frame(cbind(study, STUDY, X, y, Z))
    dat = cbind(X, y)
    X1 = X[study == 1 & Z == 0]
    y1 = y[study == 1 & Z == 0]
    dat1 = cbind(X1, y1)
    mu = apply(dat1, 2, mean)
    cov = var(dat1)
    d = apply(dat, 1, mahalanobis, center = mu, cov = cov)
    w = 1-(d-min(d))/(max(d)-min(d))
    w[Z==1] = 1
    w1 = w[study==1 & Z==0]
    w[study==1] = 1
    eps = quantile(w1, p = 0.1)
    
    X2 = X[study == 1 | Z == 0]
    y2 = y[study == 1 | Z == 0]
    Z2 = Z[study == 1 | Z == 0]
    w02 = w[study == 1 | Z == 0]
    study2 = study[study == 1 | Z == 0]
    STUDY2 = STUDY[study == 1 | Z == 0,]
    N2 = length(y2)
    data1 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w02)
    fit1 = stan('fit'=fit, 'data'=data1, warmup=1000, iter=2000, chains=1)
    fit_ss1 = extract(fit1)
    theta1 = fit_ss1$theta
    
    w0 = rep(1, N2)
    data2 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w0)
    fit2 = stan('fit'=fit, 'data'=data2, warmup=1000, iter=2000, chains=1)
    fit_ss2 = extract(fit2)
    theta2 = fit_ss2$theta
    
    w2 = w02
    w2[w2<eps] = 0
    data4 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = w2)
    fit4 = stan('fit'=fit, 'data'=data4, warmup=1000, iter=2000, chains=1)
    fit_ss4 = extract(fit4)
    theta4 = fit_ss4$theta
    
    y1 = y2[study2 == 1]
    w1 = w2[study2 == 1]
    X1 = X2[study2 == 1]
    Z1 = Z2[study2 == 1]
    N1 = length(y1)
    data3 = list(N = N1, Nw = length(w1), y = y1, x = X1, z = Z1, w = w1)
    fit3 = stan('fit'=fit, 'data'=data3, warmup=1000, iter=2000, chains=1)
    fit_ss3 = extract(fit3)
    theta3 = fit_ss3$theta
    
    
    data5 = list(N = N2, S = S, y = y2, x = X2, z = Z2,  Z = STUDY2)
    fit5 = stan('fit'=fit0, 'data'=data5, warmup=1000, iter=2000, chains=1)
    fit_ss5 = extract(fit5)
    theta5 = fit_ss5$theta
    
    
    X0 = cbind(rep(1, length(X[study != 1 & Z == 0])), X[study != 1 & Z == 0])
    y0 = y[study != 1 & Z == 0]
    Xc = cbind(rep(1, length(X[study == 1 & Z == 0])), X[study == 1 & Z == 0])
    yc = y[study == 1 & Z == 0]
    n = length(y)
    n0 = length(y0)
    a_opt = nlminb(.5, bigG, lower = 0, upper = 1, X0 = X0, Xc = Xc, y0 = y0, yc = yc)$par
    wpp = w02
    wpp[study2!=1] = a_opt
    
    data8 = list(N = N2, Nw = length(w02), y = y2, x = X2, z = Z2, w = wpp)
    #fit8 = stan(file='GaussCT_pp_wCovariates_1x.stan', data=data8, chains=0)
    fit8 = stan('fit'=fit, 'data'=data8, warmup=1000, iter=2000, chains=1)
    fit_ss8 = extract(fit8)
    theta8 = fit_ss8$theta
    
    df = data.frame(theta = c(theta1, theta2, theta3, theta4, theta5, theta8), 
                    method = c(rep('IW',1000),rep( 'FH', 1000), 
                               rep('NP', 1000), rep('TIW', 1000),
                               rep('MAP', 1000), rep('PP', 1000)))
    df_summary = df %>%
      group_by(method) %>%
      dplyr::summarize(meanTheta = mean(theta), lowTheta = quantile(theta, 0.025), 
                       uppTheta = quantile(theta, 0.975), RMSE = sqrt(mean((theta - theta_t)^2)))
    sim = rep(k, nrow(df_summary))
    df_sim = data.frame(rbind(df_sim, cbind(sim, df_summary)))
  }
  return(df_sim)
}


cases = data.frame(delta_h1 = c(0, 0, 0, 0, 5, 1),
                   delta_h2 = c(0, 0, 5, 0, 5, 1),
                   mu_h1 = c(1, 1, 1, 4, 1, 1),
                   mu_h2 = c(1, 4, 1, 4, 1, 1),
                   scn = c('exch', 'p.exch1', 'p.exch2', 'n.exch1', 'n.exch2', 'n.exch3'))
Nc = c(25, 50, 100)
cases = cbind(cases[rep(1:6, each = 3),], N = rep(Nc, 6))

 
sims = NULL                  
for (i in 1:nrow(cases)) {
  sims = rbind(sims, sim_fun(delta_h1 = cases$delta_h1[i], delta_h2 = cases$delta_h2[i], 
                             mu_h1= cases$mu_h1[i], mu_h2 = cases$mu_h2[i], Nc = cases$N[i]))
}             

sims$scn = rep(cases$scn, each = 3000)
sims$N = rep(cases$N, each = 3000)
save(sims, file = 'new_sims05.Rdata')


#plots

# df01 = sims %>%
#   filter(method %in% c('TIW', "MAP", "PP", 'NP'))
# 
# df1 = sims %>%
#   mutate(bias = (meanTheta - 0.5), IQE = (uppTheta - lowTheta), power.I = (lowTheta>0))
# 
# df2 = sims %>%
#   mutate(bias = (meanTheta - theta), IQE = (uppTheta - lowTheta))
# ggplot(sims, aes(x = method, y = RMSE, fill = method)) + 
#   geom_boxplot(size = .35, width = .5, outlier.alpha = .5, color = 'darkgrey') +
#   scale_fill_brewer(palette = 'Set2') + facet_grid(N~scn) +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'),
#         # legend.title = element_text(size = 14, face = 'bold'),
#         # legend.text = element_text(size = 14)
#         legend.position = 'none')
# ggsave()
# ggplot(df1, aes(x = method, y = bias, fill = method)) + 
#   geom_boxplot(size = .35, width = .5, outlier.alpha = 0.5, color = 'darkgrey') +
#   scale_fill_brewer(palette = 'Set2') + facet_grid(N~scn) +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'none')
# ggplot(df1, aes(x = method, y = IQE, fill = method)) + 
#   geom_boxplot(size = .35, width = .5, outlier.alpha = 0.5, color = 'darkgrey') +
#   scale_fill_brewer(palette = 'Set2') + facet_grid(scn~N) +
#   ylab('95% credible interval width') +
#   theme(axis.text = element_text(size = 14),
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         axis.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'none')
# 
# df_summary = df1 %>% group_by(method, scn, N) %>% 
#   summarize(power = mean(power.I), bias = mean(bias), RMSE = mean(RMSE), IQE = mean(IQE))
# 
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# ggplot(df_summary, aes(x = N, y = power, linetype = method, shape = method, color = method)) + 
#   geom_point() +
#   geom_line() + 
#   scale_color_manual( values = cbPalette) +
#   facet_grid(.~scn) + theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'power05.png')
# ggplot(df_summary, aes(x = N, y = RMSE, linetype = method, shape = method, color = method)) + 
#   geom_point() +
#   geom_line() + 
#   scale_color_manual( values = cbPalette) +
#   facet_grid(.~scn) + theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'RMSE05.png')
# ggplot(df_summary, aes(x = N, y = bias, linetype = method, shape = method, color = method)) + 
#   geom_point() +
#   geom_line() + 
#   scale_color_manual( values = cbPalette) +
#   facet_grid(.~scn) + theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'bias05.png')
# ggplot(df_summary, aes(x = N, y = IQE, linetype = method, shape = method, color = method)) + 
#   scale_color_manual( values = cbPalette) +
#   geom_point() +
#   geom_line() + ylab('average width of 95% CI') +
#   facet_grid(.~scn) + theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = 'bold'))
# ggsave(file = 'IQE05.png')

