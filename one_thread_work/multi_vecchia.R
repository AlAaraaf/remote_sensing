library(GPvecchia)
library(SpatialTools)
library(gstat)
library(sp)
library(MASS)
library(parallel)
library(fields)
library(nloptr)
library(LatticeKrig)

# all the means are assumed to be zero
sub_region_division = function(dat, m, layer  = 1){
  # dat -- data matrix from a coarser layer. For now suppose it's on a grid
  # m -- divide the domain into m*m  sub regions
  # layer -- indicate which layer this is
  n = nrow(dat)
  nc = ncol(dat)
  # determine sub region index
  x = dat[,1]
  y = dat[,2]
  break.x = seq(min(x), max(x), length.out = m+1)[2:m]
  break.y = seq(min(y), max(y), length.out = m+1)[2:m]
  
  inter.x = findInterval(dat[,1], break.x, rightmost.closed = TRUE)
  inter.y = findInterval(dat[,2], break.y, rightmost.closed = TRUE)
  
  dat = cbind(dat, inter.x+1 + (inter.y)*m)
  
  # determine the knot for each sub region
  # 0 means not a knot
  # j means being the knot of sub region j
  dat = cbind(dat, rep(0, n))
  for(j in 1:(m^2)){
    dj = dat[dat[, nc+1] == j,]
    xj = sort(unique(dj[,1]))
    yj = sort(unique(dj[,2]))
    med.x = xj[length(xj)/2]
    med.y = yj[length(yj)/2]
    med.index = dat[,1] == med.x & dat[,2]==med.y
    dat[med.index, nc+2] = j
  }
  return(dat)
}

sub_region_nn = function(dat.sub, m, K, option = c("mr", "cons")){
  #dat.sub -- matrix with sub region info
  n = nrow(dat.sub)
  nc = ncol(dat.sub)
  # include the neighborhood structure in the data too
  dat.sub = cbind(dat.sub, matrix(0, nrow = n, ncol = K+1))
  if(option == "mr"){
    for(i in 1:m^2){
      dat1.index = which(dat.sub[, nc-1] == i&
                           dat.sub[, nc]==0)
      dat1 = dat.sub[dat1.index,]
      ord = order_maxmin_exact(dat1[,1:2])
      #cut = min(nrow(dat1), 9)
      #ord = c(ord[1], ord[-seq(1, cut)], ord[2:cut])
      dat1 = dat1[ord,]
      dat1.index_ord = dat1.index[ord]
      NNarray = GpGp::find_ordered_nn(dat1[, 1:2], K)
      NNarray = NNarray[,-1]
      dat.sub[dat1.index_ord, nc+1] = 1:nrow(dat1)
      dat.sub[dat1.index_ord, (nc+2):(nc+K+1)] = NNarray
      }
    }else if(option =="cons"){
      dat1.index = which(dat.sub[, nc]==0)
      dat1 = dat.sub[dat1.index, ]
      ord = order_maxmin_exact(dat1[,1:2])
      #cut = min(nrow(dat1), 9)
      #ord = c(ord[1], ord[-seq(1, cut)], ord[2:cut])
      dat1 = dat1[ord,]
      dat1.index_ord = dat1.index[ord]
      NNarray = GpGp::find_ordered_nn(dat1[, 1:2], K)
      NNarray = NNarray[,-1]
      dat.sub[dat1.index_ord, nc+1] = 1:nrow(dat1)
      dat.sub[dat1.index_ord, (nc+2):(nc+K+1)] = NNarray
    }
  
  return(dat.sub)
}


sub_region_vecchia = function(dat1, knots, pars, K){
  # dat1 --sub region data matrix including conditioning set info
  # knots -- the knots one from each of the sub region
  # K-- the number of neighbors to condition on after maximin ordering
  # pars -- parameters to calculate likelihood
  
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  n = nrow(dat1)
  nc = ncol(dat1)
  M = nrow(knots)
  
  #the 6th col of dat1 records the maxmin order
  ord  = order(dat1[,6])
  dat1 = dat1[ord,]
  dat.aug  = rbind(knots, dat1)
  NNarray = dat1[,7:nc] + M
  
  log.vecchia = sapply((M+1):(M+n), function(i){
    loglik_cond(dat.aug, i, pars, M, NNarray)})
  
  sum(log.vecchia)
}



loglik_cond = function(dat.aug, pt, pars, M, NNarray){
  # dat.aug -- the knot data matrix and whole  sub region data matrix
  # pt -- the index to calculate conditional density on
  # pars -- parameters
  # M -- number of knots
  # NNarray -- matrix giving the nearest neighbors
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  stopifnot(pt > M)
  if(pt == M+1){
    cond.set = 1:M
  }else{
    cond.set = c(1:M, NNarray[pt-M, ])
    cond.set = cond.set[!is.na(cond.set)]
  }
  S11 = sigmasq + tausq
  D12 = dist2(t(dat.aug[pt, 1:2]), dat.aug[cond.set, 1:2])
  S12 = MaternFun(D12, c(sigmasq, phi, kappa))
  D22 = dist1(dat.aug[cond.set, 1:2])
  S22 = MaternFun(D22, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D22)))
  e = eigen(S22)
  S22.inv = e$vectors%*%(diag(1/e$values))%*%t(e$vectors)
  cond.mu = S12%*%S22.inv%*%dat.aug[cond.set, 3]
  cond.var  = S11 - S12%*%S22.inv%*%t(S12)
  
  -1/2*log(2*pi) -1/2*log(cond.var) - 1/2*(dat.aug[pt,3]-cond.mu)^2/cond.var
}



loglik = function(y,S){
  e = eigen(S)
  S.inv = e$vectors%*%(diag(1/e$values))%*%t(e$vectors)
  log.det = sum(log(e$values))
  -length(y)/2*log(2*pi) -1/2*log.det - 1/2*t(y)%*%S.inv%*%y
}


loglik_mr = function(pars, dat.nn, m, K){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  knots = dat.nn[dat.nn[, 5]!=0,]
  sub.loglik = 0
  for(i in 1:m^2){
    dat1 = dat.nn[dat.nn[, 4] == i&
                     dat.nn[, 5]==0,]
    sub.loglik = sub.loglik + sub_region_vecchia(dat1, knots, pars, K)
  }
  D = dist1(as.matrix(knots[, 1:2]))
  V.knots = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
  knots.loglik = loglik(y = knots[,3], V.knots)
  sub.loglik + knots.loglik
}

loglik_cons = function(pars, dat.nn, m, K){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  knots = dat.nn[dat.nn[, 5]!=0,]
  dat1 = dat.nn[dat.nn[, 5]==0,]
  
  loglik1 = sub_region_vecchia(dat1, knots, pars, K)
  
  D = dist1(as.matrix(knots[, 1:2]))
  V.knots = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
  knots.loglik = loglik(y = knots[,3], V.knots)
  
  loglik1 + knots.loglik
}

ng_loglik_mr = function(pars, dat.nn, m, K){
  -1*loglik_mr(pars, dat.nn, m, K)
}

ng_loglik_cons = function(pars, dat.nn, m, K){
  -1*loglik_cons(pars, dat.nn, m, K)
}

ng_loglik_gp = function(pars, dat, vecchia.approx){
  -1*vecchia_likelihood(dat[,3],
                        vecchia.approx, 
                        covparms = pars[1:3], 
                        nuggets = pars[4])
}

ng_loglik = function(pars, dat.nn){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  if(nrow(dat.nn)<=5000){
    D = dist1(as.matrix(dat.nn[, 1:2]))
    V = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
    -1*loglik(dat.nn[,3], V)
  }else{
    ng_loglik_cons(pars, dat.nn, m = 10, K = 40)
  }
}

sigmafn  = function(loc){
  5*(loc[1]^2+loc[2]^2+1)
}

taufn = function(loc){
  0.2*(loc[1]+1)
}


sim = function(N, m, K, pars){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  x = (1:N)/N
  y = (1:N)/N
  coords = expand.grid(x, y)
  names(coords) = c("x", "y")
  n = N*N

  obj =  circulantEmbeddingSetup(grid = list(x=x, y=y),
                                 cov.args=list(Covariance="Matern",
                                               aRange=phi,
                                               smoothness=kappa))
  z = sqrt(sigmasq)*circulantEmbedding(obj)
  #image.plot(x, y, z)
  z = as.vector(z) + sqrt(tausq)*rnorm(n)
  dat = cbind(as.matrix(coords), z)
  

  # our vecchia
  time.mr = Sys.time()
  dat.sub = sub_region_division(dat, m)
  dat.nn = sub_region_nn(dat.sub, m, K, option = "mr")
  # loglik_true_mr = loglik_mr(pars, dat.nn, m, K)
  pars_mr = optim(par = c(2, 0.04, 0.2, 0.1), fn = ng_loglik_mr,
                  dat.nn = dat.nn, m = m, K = K,
                  method = "L-BFGS-B",
                  lower = c(1, 0.001, 0.1, 0.01),
                  upper = c(20, 0.1, 10, 2),
                  control  = list(factr = 1e10))
  time.mr = Sys.time() - time.mr
  
  
  # GPvecchia pkg, conditioning on both latent and obs
  time.gp = Sys.time()
  vecchia.approx = vecchia_specify(dat[,1:2], m = K, ordering = "maxmin",
                                   cond.yz = "SGV", conditioning = "NN")
  #SGV is sparse general vecchia, condition on both obs and latent
  # loglik_true_gp = vecchia_likelihood(dat[,3], vecchia.approx,
  #                                     covparms = pars[1:3], nuggets = pars[4])
  
  # GPvecchia pkg, conditioning only on obs
  # vecchia.approx1 = vecchia_specify(dat[,1:2], m = K, ordering = "maxmin",
  #                                  cond.yz = "z", conditioning = "NN")
  # loglik_true_gp1 = vecchia_likelihood(dat[,3], vecchia.approx1,
  #                                     covparms = pars[1:3], nuggets = pars[4])
  
  pars_gp = optim(par = c(2, 0.04, 0.2, 0.1), fn = ng_loglik_gp,
                  dat = dat, vecchia.approx = vecchia.approx,
                  method = "L-BFGS-B",
                  lower = c(1, 0.001, 0.1, 0.01),
                  upper = c(20, 0.1, 10, 2),
                  control  = list(factr = 1e10))
  time.gp = Sys.time() - time.gp
  
  
  # our conservative method
  time.cons = Sys.time()
  dat.sub = sub_region_division(dat, m)
  dat.nn = sub_region_nn(dat.sub, m, K, option = "cons")
  # loglik_true_cons = loglik_cons(pars, dat.nn, m, K)
  pars_cons = optim(par = c(2, 0.04, 0.2, 0.1), fn = ng_loglik_cons,
                    dat = dat.nn, m = m, K = K,
                    method = "L-BFGS-B",
                    lower = c(1, 0.001, 0.1, 0.01),
                    upper = c(20, 0.1, 10, 2),
                    control  = list(factr = 1e10))
  time.cons = Sys.time() - time.cons
  

  #true likelihood
  time.mle = Sys.time()
  dat.sub = sub_region_division(dat, m)
  dat.nn = sub_region_nn(dat.sub, m = 10, K = 40, option = "cons")
  # loglik_true = -1*ng_loglik(pars, dat.nn)
  pars_mle = optim(par = c(2, 0.04, 0.2, 0.1), fn = ng_loglik,
                     dat.nn = dat.nn,
                     method = "L-BFGS-B",
                     lower = c(1, 0.001, 0.1, 0.01),
                     upper = c(20, 0.1, 10, 2),
                     control  = list(factr = 1e10))
  time.mle = Sys.time() - time.mle
  
  
  # time.mle = Sys.time()
  # pars_mle = nloptr(x0 = c(2, 0.04, 0.2, 0.1), 
  #                   eval_f = ng_loglik,
  #                   lb = c(1, 0.001, 0.1, 0.01),
  #                   ub = c(20, 0.1, 10, 2),
  #                   opts = list("algorithm"="NLOPT_LN_COBYLA",
  #                               "xtol_rel"=1.0e-8,
  #                               "maxeval" = 1000),
  #                     dat.nn = dat.nn
  #                     )
  # time.mle = Sys.time() - time.mle
  
  # reletive_loglik = c(loglik_true_mr/loglik_true,
  #                     loglik_true_cons/loglik_true,
  #                     loglik_true_gp/loglik_true)
  
  bias_mr = abs((pars_mr$par - pars)/pars)
  bias_gp = abs((pars_gp$par - pars)/pars)
  bias_cons = abs((pars_cons$par - pars)/pars)
  bias_mle = abs((pars_mle$par - pars)/pars)
  
  # return(list(reletive_loglik, time.mr, time.gp, time.cons, time.mle))
  
  return(list(pars_mr, pars_gp, pars_cons, pars_mle,
              bias_mr, bias_gp, bias_cons, bias_mle,
              time.mr, time.gp, time.cons, time.mle))
}

num.core = detectCores()

sigmasq = 10
kappa = 2.25
tausq = 0.25
m = 3
K = 10
set.seed(10723)
for(N in c(40, 100)){
  # when there are more data points, shrink the range s.t. it's like an increasing domain asypt
  res = sim(N, m, K, pars = c(sigmasq, 1.5/N, kappa, tausq))
  # res1 = mclapply(1:30, function(i){sim(N, m, K, pars = c(sigmasq, 1.5/N, kappa, tausq))}, mc.cores = num.core)
  # res2 = mclapply(1:30, function(i){sim(N, m, K, pars = c(sigmasq, 3/N, kappa, tausq))}, mc.cores = num.core)
  save.image(paste0("v", N, ".RData"))
  
}

res.plot = function(res, N, phi){
  res = matrix(unlist(res), ncol = 20, byrow = TRUE)
  boxplot.matrix(res[,c(1,5,9,13)], main = paste("percent bias on psill",
                                               "N = ", N,
                                               "range = ", phi, sep = " "))
  boxplot.matrix(res[,c(2,6,10,14)], main = paste("percent bias on range",
                                                "N = ", N,
                                                "range = ", phi, sep = " "))
  boxplot.matrix(res[,c(3,7,11,15)], main = paste("percent bias on smooth",
                                                "N = ", N,
                                                "range = ", phi, sep = " "))
  boxplot.matrix(res[,c(4,8,12,16)], main = paste("percent bias on nugget",
                                                "N = ", N,
                                                "range = ", phi, sep = " "))
  boxplot.matrix(res[,17:20], main = paste("time",
                                            "N = ", N,
                                            "range = ", phi, sep = " "))
}
N=40
res.plot(res1, N, 0.4/N)
res.plot(res2, N, 1/N)


