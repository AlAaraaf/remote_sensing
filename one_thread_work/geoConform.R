library(SpatialTools)
library(gstat)
library(sp)
library(MASS)
library(lava)
library(parallel)
library(tidyr)

upper.quantile = function(v, prob, w=NULL, sorted=FALSE) {
  if (is.null(w)) w = rep(1,length(v))
  if (!sorted) { o = order(v); v = v[o]; w = w[o] }
  i = which(cumsum(w/sum(w)) >= prob)
  if (length(i)==0) return(Inf)
  else return(v[min(i)])
}

loess.cv = function(x, y, K = 5, rep = 1, span.range = seq(0.05, 1, by = 0.05)){
  d = data.frame(x, y)
  colnames(d)[-1] = 'y'
  K = as.integer(K)
  if(K < 1){
    stop("Fold number must be a positive integer")
  }
  else{
    folds = foldr(nrow(d), K, rep)
  }
  cv.mae = lapply(span.range, function(s) measure.cv(s, d, K, folds))
  best.span = span.range[which.min(cv.mae)]
  mod.best = loess(d$y ~ ., data = d, degree = 1, span = best.span, 
                   control = loess.control(surface = "direct"))
  return(mod.best)
}

measure.cv = function(s, data, K, folds){
  mabse = 0
  for(k in 1:K){
    fold = as.vector(folds[[1]][[k]])
    dtest = data[fold, ]
    dtrain  = data[-fold, ]
    model.k  = loess(dtrain$y ~ ., data = dtrain, degree = 1, span = s, 
                     control = loess.control(surface = "direct"))
    yhat.k  = predict(model.k, dtest)
    mabse = mabse + sum(abs(yhat.k - dtest$y))
  }
  mabse = mabse/nrow(data)
  return(mabse)
}

# dist.mx = function(coords1, coords2){
#   coords1 = as.matrix(coords1)
#   coords2 = as.matrix(coords2)
#   n1 = nrow(coords1)
#   n2 = nrow(coords2)
#   sapply(1:n1,function(i){
#     sapply(1:n2,function(j){
#       dst(coords1[i,],coords2[j,])
#     })
#   })
# }
# 
# dst = function(s1, s2){
#   sqrt(sum((s1-s2)^2))
# }

loglik.fix = function(par,data, sigma.e, d.mtx = NULL){
  #sigma.e is a fitted curve evaluated at coords, should input sigma^2(s)
  sigmasq = par[1]
  phi = par[2]
  mu = par[3]
  
  coords = as.matrix(data[c("x", "y")])
  y = as.matrix(data["data"])
  n = nrow(coords)
  
  if(is.null(d.mtx)){
    d.mtx = dist1(coords)
  }
  
  S = sigmasq*exp(-d.mtx/phi) + diag(sigma.e)
  S.inv = solve(S)
  e = eigen(S)
  log.det = sum(log(e$values))
  l = -1/2*log.det - 1/2*t(y-mu)%*%S.inv%*%(y-mu)
  
  return(-1*l)
}

grd.fix = function(par,data, sigma.e, d.mtx = NULL){
  #sigma.e is a fitted curve evaluated at coords
  sigmasq = par[1]
  phi = par[2]
  mu = par[3]
  
  coords = as.matrix(data[c("x", "y")])
  y = as.matrix(data["data"])
  n = nrow(coords)
  
  if(is.null(d.mtx)){
    d.mtx = dist1(coords)
  }
  
  S = sigmasq*exp(-d.mtx/phi) + diag(sigma.e)
  S.inv = solve(S)
  
  dldm = sum(S.inv%*%(y-mu))
  
  # dldS = S.inv +1/2*S.inv*diag(n) +S.inv%*%(y-mu)%*%t(y-mu)%*%S.inv-
  #   1/2*S.inv%*%(y-mu)%*%t(y-mu)%*%S.inv*diag(n)
  dldS = -1/2*S.inv + 1/2*S.inv%*%(y-mu)%*%t(y-mu)%*%S.inv
  dSds = exp(-d.mtx/phi)
  dSdp = sigmasq*exp(-d.mtx/phi)*d.mtx/phi^2
  
  dlds = sum(dldS*dSds)
  dldp = sum(dldS*dSdp)
  
  return(-1*c(dlds, dldp, dldm))
  
}

## kriging

## local kriging
#library(laGP)

## GSCP
gscp.conform = function(train.d, newlocn){
  n = nrow(train.d)
  n0 = nrow(newlocn)
  samp = sample(1:n, size = n/2)
  prop.d = train.d[samp, ]
  cali.d = train.d[-samp, ]
  

  pred.loc = rbind(cali.d[c("x","y")], newlocn)
 
  D = dist1(as.matrix(prop.d[c("x", "y")]))
  D.x = dist2(as.matrix(prop.d[c("x", "y")]), as.matrix(pred.loc))
  C.x = coincident(as.matrix(prop.d[c("x", "y")]), as.matrix(pred.loc))
  D.p = dist1(as.matrix(pred.loc))
  
  coordinates(prop.d) = ~x+y
  vario = gstat::variogram(data~1, prop.d)
  vario.fit = fit.variogram(vario,model = vgm("Exp"))
  sigmasq = vario.fit$psill[2]
  tausq = vario.fit$psill[1]
  phi = vario.fit$range[2]
  
  V = sigmasq*exp(-D/phi)+diag(rep(tausq,n/2))
  V.x = sigmasq*exp(-D.x/phi)
  if(nrow(C.x)!=0){
    for(i in 1: nrow(C.x)){
      c.x = C.x[i,]
      V.x[c.x[1], c.x[2]] = V.x[c.x[1], c.x[2]]+ tausq
    }
  }
  
  V.p = sigmasq*exp(-D.p/phi)+diag(rep(tausq,n/2+n0))
  k = krige.ok(prop.d$data, V,V.p, V.x)
  
  r = abs(k$pred- c(cali.d$data, rep(0, n0)))
  scores = r/(k$mspe+1e-15)
  cali.scores = scores[1:(n/2)]
  test.pred = k$pred[(n/2+1):(n/2+n0)]
  test.var = k$mspe[(n/2+1):(n/2+n0)]+1e-15
  return(list(scores = cali.scores, test.pred = test.pred, test.var = test.var))
}

gscp.pi = function(alpha, scores, pred, v){
  q = upper.quantile(c(scores, Inf), 1-alpha)
  lower = pred-q*v
  upper = pred+q*v
  PI = cbind(lower, upper)
  return(PI)
}

## LSCP
lscp.conform = function(train.d, newlocn, eta = 2){
  n = nrow(train.d)
  n0 = nrow(newlocn)
  samp = sample(1:n, size = n/2)
  prop.d = train.d[samp, ]
  cali.d = train.d[-samp, ]
  
  pred.loc = rbind(cali.d[c("x","y")], newlocn)
  
  D = dist1(as.matrix(prop.d[c("x", "y")]))
  D.x = dist2(as.matrix(prop.d[c("x", "y")]), as.matrix(pred.loc))
  C.x = coincident(as.matrix(prop.d[c("x", "y")]), as.matrix(pred.loc))
  D.p = dist1(as.matrix(pred.loc))
  
  coordinates(prop.d) = ~x+y
  vario = gstat::variogram(data~1, prop.d)
  vario.fit = fit.variogram(vario, model = vgm("Exp"))
  sigmasq = vario.fit$psill[2]
  tausq = vario.fit$psill[1]
  phi = vario.fit$range[2]
  
  V = sigmasq*exp(-D/phi)+diag(rep(tausq,n/2))
  V.x = sigmasq*exp(-D.x/phi)
  if(nrow(C.x)!=0){
    for(i in 1: nrow(C.x)){
      c.x = C.x[i,]
      V.x[c.x[1], c.x[2]] = V.x[c.x[1], c.x[2]]+ tausq
    }
  }
  
  V.p = sigmasq*exp(-D.p/phi)+diag(rep(tausq,n/2+n0))
  
  rad = 2*eta
  M = 15
  nn = c()
  for(i in 1: ncol(D.x)){
    m = sum(D.x[,i]<=rad)
    nn = c(nn, m)
  }
  M = max(M, min(nn))
  
  k = data.frame()
  for(i in 1:(n/2+n0)){
    d.x = D.x[,i]
    o = order(d.x)[1:M]
    y = prop.d$data[o]
    k.i = krige.ok(y, V[o,o], as.matrix(V.p[i,i]), as.matrix(V.x[o,i], ncol = 1))
    k = rbind(k, k.i)
  }

  r = abs(k$pred- c(cali.d$data, rep(0, n0)))
  scores = r/(k$mspe+1e-15)
  cali.scores = scores[1:(n/2)]
  test.pred = k$pred[(n/2+1):(n/2+n0)]
  test.var = k$mspe[(n/2+1):(n/2+n0)]+1e-15
  return(list(scores = cali.scores, 
              test.pred = test.pred, 
              test.var = test.var,
              num = M,
              cali.d = cali.d))
}

lscp.pi = function(alpha, scores, pred, v,cali.d, newlocn, M){
  cali.coords = cali.d[c("x", "y")]
  D.x = dist2(as.matrix(cali.coords), as.matrix(newlocn))
  PI = data.frame()
  M = max(M, ceiling(2/alpha))
  M = min(M, nrow(cali.d))
  for(i in 1:nrow(newlocn)){
    d.x = D.x[,i]
    o = order(d.x)[1:M]
    q = upper.quantile(c(scores[o], Inf), 1-alpha)
    lower = pred[i]-q*v[i]
    upper = pred[i]+q*v[i]
    PI = rbind(PI, cbind(lower, upper))
  }
  return(PI)
}

## sLSCP

slscp.pi = function(alpha, scores, pred, v, cali.d, newlocn, M, eta = 2){
  cali.coords = cali.d[c("x", "y")]
  D.x = dist2(as.matrix(cali.coords), as.matrix(newlocn))
  PI = data.frame()
  M = max(M, ceiling(2/alpha))
  M = min(M, nrow(cali.d))
  for(i in 1:nrow(newlocn)){
    d.x = D.x[,i]
    o = order(d.x)[1:M]
    w = exp(-d.x[o]^2/2/eta^2)
    q = upper.quantile(c(scores[o], Inf), 1-alpha, c(w, 1))
    lower = pred[i]-q*v[i]
    upper = pred[i]+q*v[i]
    PI = rbind(PI, cbind(lower, upper))
  }
  return(PI)
}


## ASCP
ascp.conform = function(train.d, newlocn){
  n = nrow(train.d)
  n0 = nrow(newlocn)
  samp = sample(1:n, size = n/2)
  prop.d = train.d[samp, ]
  cali.d = train.d[-samp, ]
  
  coords.train = prop.d[c("x", "y")]
  D = dist1(as.matrix(coords.train))
  #estimate sigma.e by differencing (piecewise constant + smoothing)
  k = 4 #use 4 nearest neighbors to calculate the difference sequence
  l = 20 
  ######
  #can add a threshold for max distance
  #######
  const = matrix(0, nrow = n/2, ncol = l)
  diff = rep(0, n/2)
  sigma.d = rep(0, n/2)
  for(i in 1:(n/2)){
    d = D[i,]
    o = order(d)
    const[i,] = o[2:(l+1)]
    nei = o[2:(k+1)]
    diff[i] = (4*prop.d$data[i] - sum(prop.d$data[nei]))^2/20
  }
  for(i in 1:n/2){
    sigma.d[i] = mean(diff[c(i,const[i,])])
  }
  
 
  #first optimization
  #determine initial values
  m = mean(prop.d$data)
  prop.g = prop.d
  coordinates(prop.g) = ~x+y
  v.sample = gstat::variogram(data~1, prop.g)
  v = mean(tail(v.sample)$gamma)
  r0 = max(v.sample$dist)/3
  theta0 = optim(par = c(v, r0, m), fn = loglik.fix, gr = grd.fix, 
                 data = prop.d, sigma.e = sigma.d, d.mtx = D, 
                 method = "L-BFGS-B", 
                 lower = c(v/5, r0/5, -10*abs(m)), 
                 upper = c(v*100, r0*100, 10*abs(m)),
                 control = list(trace = 6))
  
  #kriging with measurement error for the signal, on the proper training set
  sigmasq = theta0$par[1]
  phi = theta0$par[2]
  mu = theta0$par[3]
  
  S.z = sigmasq*exp(-D/phi)
  S = S.z + sigma.d*diag(n/2)
  # S.inv = solve(S)
  # lambda = S.inv%*%S.z + S.inv%*%rep(1, n/2)%*%solve(t(rep(1, n/2))%*%S.inv%*%rep(1, n/2))%*%(rep(1, n/2)-t(rep(1, n/2))%*%S.inv%*%S.z)
  # pred.z = t(lambda)%*%prop.d$data
  pred.z = rep(0, n/2)
  for(i in 1:(n/2)){
    tau.z = S.z[i, ]
    tau.z = tau.z[-i]
    S.i = S[-i, -i]
    S.inv.i  = solve(S.i)
    lambda = S.inv.i%*%tau.z + S.inv.i%*%rep(1, n/2-1)%*%solve(t(rep(1, n/2-1))%*%S.inv.i%*%rep(1, n/2-1))%*%(1-t(rep(1, n/2-1))%*%S.inv.i%*%tau.z)
    pred.z[i] = t(lambda)%*%prop.d$data[-i]
  }
  
  #fit s(x,y)
  r = abs(pred.z-prop.d$data)
  s_mod = smooth.spline(prop.d$x, r)
  sigma.e = predict(s_mod,prop.d$x)$y
  # s_mod = loess.cv(prop.d$x, r)
  # sigma.e = predict(s_mod, prop.d$x)
  
  # maximize the log-likelihood with fitted sigma.e
  theta1 = optim(par = c(sigmasq, phi, mu), fn = loglik.fix,gr = grd.fix,
                 data = prop.d, sigma.e = sigma.e^2, d.mtx = D,
                 method = "L-BFGS-B", 
                 lower = c(v/5, r0/5, -10*abs(m)), 
                 upper = c(v*100, r0*100, 10*abs(m)), 
                 control = list(trace = 6))
  sigmasq = theta1$par[1]
  phi = theta1$par[2]
  mu = theta1$par[3]
  
  #calculate nonconformity scores
  S = sigmasq*exp(-D/phi) + diag(sigma.e^2)
  S.inv = solve(S)
  D.x = dist2(as.matrix(coords.train), as.matrix(rbind(cali.d[c("x","y")], newlocn)))
  tau.z = sigmasq*exp(-D.x/phi)
  lambda = S.inv%*%tau.z + S.inv%*%rep(1, n/2)%*%solve(t(rep(1, n/2))%*%S.inv%*%rep(1, n/2))%*%(rep(1, n/2+n0)-t(rep(1, n/2))%*%S.inv%*%tau.z)
  pred.z = t(lambda)%*%prop.d$data
  
  r = abs(as.vector(pred.z) - c(cali.d$data, rep(0, n0)))
  sigma.e0 = predict(s_mod, c(cali.d$x, newlocn$x))$y
  
  scores = r/sigma.e0
  cali.scores = scores[1:(n/2)]
  test.pred = as.vector(pred.z)[(n/2+1):(n/2+n0)]
  test.sigma.e = sigma.e0[(n/2+1):(n/2+n0)]
  return(list(scores = cali.scores, test.pred = test.pred, test.sigma.e = test.sigma.e,
              train.x = prop.d$x, differencing = sigma.d, smoothing = sigma.e^2,
              theta0 = theta0$par, theta1 = theta1$par))
}

ascp.pi = function(alpha, scores, pred, sigma.e){
  q = upper.quantile(c(scores, Inf), 1-alpha)
  lower = pred-q*sigma.e
  upper = pred+q*sigma.e
  PI = cbind(lower, upper)
  return(PI)
}

######### data generation
s1 = function(x, y){
  return((x+3)/10)
}

s2 = function(x,y){
  return(0.3+exp(-(x-5)^2))
}

sigmasq = 4
phi = 10
mu = 0
N = 30
N0 = 10

sim = function(alpha=0.05, sigmasq = 2, phi = 10, mu = 0, N = 30, N0 = 10, s){
  x = c(1:N)/(N/10)
  y = c(1:N)/(N/10)
  coords = expand.grid(x, y)
  names(coords) = c("x", "y")
  
  x0 = runif(N0, 0, 10)
  y0 = runif(N0, 0, 10)
  newlocn = expand.grid(x0, y0)
  names(newlocn) = c("x", "y")
  
  loc  = rbind(coords, newlocn)
  # coordinates(loc) = ~x+y
  # geo.dummy = gstat(formula = z~1, dummy = TRUE, beta = mu, model = vgm(psill = sigmasq, "Exp", range = phi, nugget = 0), nmax = 200)
  # geo.z = predict(geo.dummy, loc, nsim = 1)
  n = N*N
  n0 = N0*N0
  D.full = dist1(as.matrix(loc))
  V.full = sigmasq*exp(-D.full/phi)
  geo.z = rmvnorm(nsim = 1, mu = rep(mu, n+n0), V = V.full)
  
  train.z = geo.z[1:n]
  test.z = geo.z[(n+1):(n+n0)]
  
  train.e = rnorm(n, 0, sd = 1)
  test.e = rnorm(n0, 0, sd = 1)
  
  train.y = mu + train.z + s(coords$x,coords$y)*train.e
  test.y = mu + test.z + s(newlocn$x, newlocn$y)*test.e
  
  train.d = data.frame(coords, train.y) 
  colnames(train.d) = c("x", "y", "data")
  
  gscp = gscp.conform(train.d, newlocn)
  lscp = lscp.conform(train.d, newlocn)
  ascp = ascp.conform(train.d, newlocn)
  
  cali.d = lscp$cali.d
  
  g.pi = gscp.pi(alpha, gscp$scores, gscp$test.pred, gscp$test.var)
  l.pi = lscp.pi(alpha, lscp$scores, lscp$test.pred, lscp$test.var, 
                 cali.d, newlocn, lscp$num)
  sl.pi = slscp.pi(alpha, lscp$scores, lscp$test.pred, lscp$test.var, 
                 cali.d, newlocn, lscp$num)
  a.pi = ascp.pi(alpha, ascp$scores, ascp$test.pred, ascp$test.sigma.e)
  res = data.frame(x = newlocn$x,
                   g.cov = g.pi[,1] <=test.y & test.y <= g.pi[,2], 
                   g.len = g.pi[,2]-g.pi[,1],
                   l.cov = l.pi[,1] <=test.y & test.y <= l.pi[,2], 
                   l.len = l.pi[,2]-l.pi[,1],
                   sl.cov = sl.pi[,1] <=test.y & test.y <= sl.pi[,2], 
                   sl.len = sl.pi[,2]-sl.pi[,1],
                   a.cov = a.pi[,1] <=test.y & test.y <= a.pi[,2], 
                   a.len = a.pi[,2]-a.pi[,1])
  res.inter = list(x = ascp$train.x, 
                   differencing  = ascp$differencing,
                   smoothing = ascp$smoothing,
                   theta0 = ascp$theta0,
                   theta1 = ascp$theta1)
  return(list(res = res, res.inter = res.inter))
}

sim.ngp = function(alpha=0.05, sigmasq = 2, phi = 10, mu = 0, N = 30, N0 = 10, s){
  x = c(1:N)/(N/10)
  y = c(1:N)/(N/10)
  coords = expand.grid(x, y)
  names(coords) = c("x", "y")
  
  x0 = runif(N0, 0, 10)
  y0 = runif(N0, 0, 10)
  newlocn = expand.grid(x0, y0)
  names(newlocn) = c("x", "y")
  
  loc  = rbind(coords, newlocn)
  # coordinates(loc) = ~x+y
  # geo.dummy = gstat(formula = z~1, dummy = TRUE, beta = mu, model = vgm(psill = sigmasq, "Exp", range = phi, nugget = 0), nmax = 200)
  # geo.z = predict(geo.dummy, loc, nsim = 1)
  n = N*N
  n0 = N0*N0
  D.full = dist1(as.matrix(loc))
  V.full = sigmasq*exp(-D.full/phi)
  geo.z = rmvnorm(nsim = 1, mu = rep(mu, n+n0), V = V.full)
  
  train.z = geo.z[1:n]
  test.z = geo.z[(n+1):(n+n0)]
  
  train.e = rnorm(n, 0, sd = 1)
  test.e = rnorm(n0, 0, sd = 1)
  
  train.y = mu + (train.z)^2 + s(coords$x,coords$y)*train.e
  test.y = mu + (test.z)^2 + s(newlocn$x, newlocn$y)*test.e
  
  train.d = data.frame(coords, train.y) 
  colnames(train.d) = c("x", "y", "data")
  
  gscp = gscp.conform(train.d, newlocn)
  lscp = lscp.conform(train.d, newlocn)
  ascp = ascp.conform(train.d, newlocn)
  
  cali.d = lscp$cali.d
  
  g.pi = gscp.pi(alpha, gscp$scores, gscp$test.pred, gscp$test.var)
  l.pi = lscp.pi(alpha, lscp$scores, lscp$test.pred, lscp$test.var, 
                 cali.d, newlocn, lscp$num)
  sl.pi = slscp.pi(alpha, lscp$scores, lscp$test.pred, lscp$test.var, 
                   cali.d, newlocn, lscp$num)
  a.pi = ascp.pi(alpha, ascp$scores, ascp$test.pred, ascp$test.sigma.e)
  res = data.frame(x = newlocn$x,
                   g.cov = g.pi[,1] <=test.y & test.y <= g.pi[,2], 
                   g.len = g.pi[,2]-g.pi[,1],
                   l.cov = l.pi[,1] <=test.y & test.y <= l.pi[,2], 
                   l.len = l.pi[,2]-l.pi[,1],
                   sl.cov = sl.pi[,1] <=test.y & test.y <= sl.pi[,2], 
                   sl.len = sl.pi[,2]-sl.pi[,1],
                   a.cov = a.pi[,1] <=test.y & test.y <= a.pi[,2], 
                   a.len = a.pi[,2]-a.pi[,1])
  res.inter = list(x = ascp$train.x, 
                   differencing  = ascp$differencing,
                   smoothing = ascp$smoothing,
                   theta0 = ascp$theta0,
                   theta1 = ascp$theta1)
  return(list(res = res, res.inter = res.inter))
}


num.core = detectCores()

res1 = mclapply(1:100, function(x){sim(s = s1)}, mc.cores = num.core)
res2 = mclapply(1:100, function(x){sim(s = s2)}, mc.cores = num.core)
res3 = mclapply(1:100, function(x){sim.ngp(s = s1)}, mc.cores = num.core)
res4 = mclapply(1:100, function(x){sim.ngp(s = s2)}, mc.cores = num.core)
res5 = mclapply(1:100, function(x){sim(alpha=0.05, sigmasq = 0.5, phi = 10,
                                       mu = 0, N = 30, N0 = 10, s=s1)}, mc.cores = num.core)
res6 = mclapply(1:100, function(x){sim(alpha=0.05, sigmasq = 0.5, phi = 5,
                                       mu = 0, N = 30, N0 = 10, s=s1)}, mc.cores = num.core)
save.image("geoConform.RData")


extract.result = function(res, alpha=0.05){
  #result plots
  experiment_bool = data.frame()
  len = data.frame()
  cov = data.frame()
  #diagnose the intermediate steps
  inter_sigma = data.frame()
  theta0 = data.frame()
  theta1 = data.frame()
  for(i in  1:length(res)){
    exp_res = res[[i]]$res
    experiment_bool = rbind(experiment_bool, exp_res[c("x", 
                                                       "g.cov", 
                                                       "l.cov", 
                                                       "sl.cov", 
                                                       "a.cov")])
    len = rbind(len, cbind(mean(exp_res$g.len), 
                           mean(exp_res$l.len), 
                           mean(exp_res$sl.len), 
                           mean(exp_res$a.len)))
    cov = rbind(cov, cbind(mean(exp_res$g.cov), 
                           mean(exp_res$l.cov), 
                           mean(exp_res$sl.cov), 
                           mean(exp_res$a.cov)))
    
    inter_res = res[[i]]$res.inter
    inter_sigma = rbind(inter_sigma, cbind(inter_res$x,
                                           inter_res$differencing,
                                           inter_res$smoothing) )
    theta0 = rbind(theta0, inter_res$theta0)
    theta1 = rbind(theta1, inter_res$theta1)
  }
  colnames(cov) = c("GSCP", "LSCP", "sLSCP", "ASCP")
  colnames(len) = c("GSCP", "LSCP", "sLSCP", "ASCP")
  colnames(inter_sigma) = c("x", "diff", "smooth")
  colnames(theta0) = c("sigmasq", "phi", "mu")
  colnames(theta1) = c("sigmasq", "phi", "mu")
  theta0["ratio"] = theta0$sigmasq/theta0$phi
  theta1["ratio"] = theta1$sigmasq/theta1$phi
  
  par(mar=c(5.1, 4.1, 4.1, 2.1), xpd = FALSE)
  boxplot(cov, main = "marginal coverage")
  lines(c(0:5),rep(1-alpha, 6), col = "red")
  boxplot(len, main = "average interval length")
  g.line = loess.cv(experiment_bool$x, experiment_bool$g.cov)
  l.line = loess.cv(experiment_bool$x, experiment_bool$l.cov)
  sl.line = loess.cv(experiment_bool$x, experiment_bool$sl.cov)
  a.line = loess.cv(experiment_bool$x, experiment_bool$a.cov)
  par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
  x = seq(0, 10, by = 0.05)
  plot(x, predict(g.line, x), type = "l", main = "Local coverage", xlab = "x coordinate", ylab = "coverage", col = "yellow")
  lines(x, predict(l.line, x), col = "blue")
  lines(x, predict(sl.line, x), col = "green")
  lines(x, predict(a.line, x), col = "red")
  lines(x, rep(1-alpha, length(x)), col = "black", lty = 3)
  legend("topright",inset = c(-0.6,0), lty=1, col=c('black', 'yellow', 'blue', 'green', 'red'),
         legend=c("Nominal","GSCP", "LSCP", "sLSCP","ASCP"))
  
  # hist(theta0$ratio, main = "theta0, sigmasq/phi")
  # points(x=0.2, y=0, pch = 20, col ="red")
  # hist(theta0$mu, main = "theta0, mu")
  # points(x=0, y=0, pch = 20, col ="red")
  # 
  # hist(theta1$ratio, main = "theta1, sigmasq/phi")
  # points(x=0.2, y=0, pch = 20, col ="red")
  # hist(theta1$mu, main = "theta1, mu")
  # points(x=0, y=0, pch = 20, col ="red")
  # x = sort(unique(inter_sigma$x))
  # for(xx in x){
  #   inter_s = inter_sigma[inter_sigma$x == xx, ]
  #   hist(inter_s$diff, main = paste("Diff res at x=", xx))
  #   points(x=s2(xx, 0)^2, y=0, pch = 20, col ="red")
  #   hist(inter_s$smooth, main = paste("Smooth res at x=", xx))
  #   points(x=s2(xx, 0)^2, y=0, pch = 20, col ="red")
  # }
  
}
extract.result(res1)
extract.result(res2)
extract.result(res3)
extract.result(res4)
extract.result(res5)
extract.result(res6)


y = rnorm(500,0,1)
s = abs(y-mean(y))
hist(s, xlab = "scores", main = "Histogram of scores", xlim = c(0,4))
text(2, 190, expression(delta[i] ==abs(Y[i]-bar(Y))))
points(x = 0.75, y=0, pch = 4, col = "blue", lwd = 5)
points(x = 3.75, y= 0, pch = 4, col ="red", lwd = 5)
text(1, 60,  "Less than 62% of scores", col = "blue")
text(3, 60, "Less than 4% of scores", col = "red")

for(i in 1:length(res6)){
  res = res6[[i]]
  if(min(res$res.inter$smoothing)<0.09){
    print(i)
  }
}

res = res6[[98]]
View(res$res[res$res$a.cov==0,])

res = res1[[71]]
View(res$res[res$res$a.cov==0,])

res = res1[[86]]
View(res$res[res$res$a.cov==0,])

X.false = c()
for(i in 1:length(res1)){
  res = res1[[i]]$res
  x = res$x[res$a.cov==0]
  X.false = c(X.false,x)
}
