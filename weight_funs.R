
# ' Multinomial/Gaussian model-based weights
#'
#' Returns weights computed based on the joint distribution of discrete and continuous variables 
#' 
#'
#' @param cat.v matrix of categorical variables 
#' @param cont.v matrix of continuous variables
#' @param s.ind study index; 1 for current, 0 for historical
#' @return vector of weights
#'
#' @export

fit_simmodel = function(cat.v, cont.v, s.ind) {
  cont.v1 = cont.v[s.ind,]
  z = apply(cat.v, 1, paste, collapse = '')
  z = as.numeric(as.factor(z))
  z1 = z[s.ind]
  lev = sort(unique(z1))
  if(is.null(dim(cont.v))) {
    mu = mean(cont.v1[z1==1])
    Sigma = var(cont.v1[z1==1])
  } else {
    if (!is.null(dim(cont.v1[z1==1,]))) {
      mu = apply(cont.v1[z1==1,], 2, mean)
      Sigma = cov(cont.v1[z1==1,])
    } else {
      mu = cont.v1[z1==1,]
      Sigma =  matrix(NA, ncol(cont.v1), ncol(cont.v1))
    }
  }
  for (i in 2:length(lev)) {
    if(is.null(dim(cont.v1))) {
      mu = c(mu, mean(cont.v1[z1==i]))
      Sigma = c(Sigma, var(cont.v1[z1==i]))
      } else {
        if (!is.null(dim(cont.v1[z1==i,]))) {
          mu = rbind(mu, apply(cont.v1[z1==i,], 2, mean))
          Sigma = Matrix::bdiag(Sigma, cov(cont.v1[z1==i,]))
        } else {
          mu = rbind(mu, cont.v1[z1==i,])
          Sigma = Matrix::bdiag(Sigma, matrix(NA, ncol(cont.v1), ncol(cont.v1)))
        }
  }}
  px = c()
  for (j in lev) {
    px[j] = mean(z1==j)
  }
  px[is.na(px)] = 0
  pars = list(mu = mu, Sigma = Sigma, px = px)
  w = rep(NA, length(z))
  for (i in 1:length(lev)) {
    if(is.null(dim(cont.v))) {
      w[z==i] = px[i]*dnorm(cont.v[z==i], mu[i], sqrt(Sigma[i]))
    } else {
      w[z==i] = px[i]*apply(cont.v[z==i, ], 1, dmvnorm, mean = mu[i,], sigma = matrix(Sigma[(2*(i-1)+1):(2*i), (2*(i-1)+1):(2*i)], 2, 2))
    }}
  w[is.na(w)] = 0  
  w[w == Inf] = 1
  w[w>1] = 1
  u = max(w[w<1])
  w[w==1] = u
  w = (w - min(w))/(max(w) - min(w))
  out = list(pars = pars, weights = w)
  return(out)
}








