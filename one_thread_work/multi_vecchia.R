library(GPvecchia)
library(SpatialTools)
library(gstat)
library(sp)
library(MASS)
library(parallel)
library(fields)
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


sub_region_vecchia = function(dat1, knots, pars, K){
  # dat1 --sub region data matrix
  # knots -- the knots one from each of the sub region
  # K-- the number of neighbors to condition on after maximin ordering
  # pars -- parameters to calculate likelihood
  
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  n = nrow(dat1)
  M = nrow(knots)
  
  ord = order_maxmin_exact(dat1[,1:2])
  dat1 = dat1[ord,]
  dat.aug  = rbind(knots, dat1)
  
  log.vecchia = sapply((M+1):(M+n), function(i){
    loglik_cond(dat.aug, i, pars, M, K)})
  
  sum(log.vecchia)
}



loglik_cond = function(dat.aug, pt, pars,M, K){
  # dat.aug -- the knot data matrix and whole  sub region data matrix
  # pt -- the index to calculate conditional density on
  # pars -- parameters
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  stopifnot(pt > M)
  if(pt == M+1){
    cond.set = 1:M
  }else{
    cond.set = c(1:M, max(M+1, pt-K):(pt-1))
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


loglik_mr = function(pars, dat, m, K){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  dat.sub = sub_region_division(dat, m)
  knots = dat.sub[dat.sub[, ncol(dat.sub)]!=0,]
  sub.loglik = 0
  for(i in 1:m^2){
    dat1 = dat.sub[dat.sub[, ncol(dat.sub)-1] == i&
                     dat.sub[, ncol(dat.sub)]==0,]
    sub.loglik = sub.loglik + sub_region_vecchia(dat1, knots, pars, K)
  }
  D = dist1(as.matrix(knots[, 1:2]))
  V.knots = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
  knots.loglik = loglik(y = knots[,3], V.knots)
  sub.loglik + knots.loglik
}

loglik_cons = function(pars, dat, m, K){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  
  dat.sub = sub_region_division(dat, m)
  knots = dat.sub[dat.sub[, ncol(dat.sub)]!=0,]
  dat1 = dat.sub[dat.sub[, ncol(dat.sub)]==0,]
  
  loglik1 = sub_region_vecchia(dat1, knots, pars, K)
  
  D = dist1(as.matrix(knots[, 1:2]))
  V.knots = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
  knots.loglik = loglik(y = knots[,3], V.knots)
  
  loglik1 + knots.loglik
}

ng_loglik_mr = function(pars, dat, m, K){
  -1*loglik_mr(pars, dat, m, K)
}

ng_loglik_cons = function(pars, dat, m, K){
  -1*loglik_cons(pars, dat, m, K)
}

ng_loglik_gp = function(pars, dat, vecchia.approx){
  -1*vecchia_likelihood(dat[,3],
                        vecchia.approx, 
                        covparms = pars[1:3], 
                        nuggets = pars[4])
}

ng_loglik = function(pars, dat){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  if(nrow(dat)<10000){
    D = dist1(as.matrix(dat[, 1:2]))
    V = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
    -1*loglik(dat[,3], V)
  }else{
    ng_loglik_cons(pars, dat, m = 10, K = 40)
  }
}

krigeSimCE1 = function (formula, data, newdata, model, n = 1, ext = 2) 
{
  stopifnot(is(model, "variogramModel"))
  stopifnot(gridded(newdata))
  if (!missing(data)) 
    stopifnot(identical(data@proj4string@projargs, newdata@proj4string@projargs))
  varName <- all.vars(formula[[2]])
  condSim <- TRUE
  if (missing(data)) {
    condSim <- FALSE
    message("[No data provided: performing unconditional simulation.]")
  }
  else {
    message("[Performing conditional simulation.]")
  }
  covMat <- gstat:::ceWrapOnTorusCalcCovRow1(newdata, model, ext = ext)
  sims <- gstat:::ceSim(covMat, n, newdata@grid@cells.dim, newdata@grid.index)
  colnames(sims) <- paste0(varName, ".sim", 1:n)
  if (!condSim) {
    if ("data" %in% slotNames(newdata)) 
      newdata@data <- cbind(newdata@data, sims)
    else newdata = addAttrToGeom(newdata, as.data.frame(sims))
    return(newdata)
  }
  obsMeanField <- krige(formula, data, newdata, model)
  simMeanObsLoc <- krigeMultiple(as.formula(paste0("var1.pred ~", 
                                                   formula[[3]])), obsMeanField, data, model, sims)
  simMeanFields <- krigeMultiple(as.formula(paste0(varName, 
                                                   "~", formula[[3]])), data, newdata, model, simMeanObsLoc)
  sims <- obsMeanField@data$var1.pred + sims - simMeanFields
  if ("data" %in% slotNames(newdata)) {
    newdata@data <- cbind(newdata@data, sims)
    return(newdata)
  }
  newdata = addAttrToGeom(newdata, as.data.frame(sims))
}

st_generator = function(loc, pars, method = c("seq", "CE"), tr = 1){
  sigmasq = pars[1]
  phi = pars[2]
  kappa = pars[3]
  tausq = pars[4]
  modVgm = vgm(sigmasq, "Mat", phi, tausq, kappa = kappa)
  if(method == "CE"){
    gridded(loc) = ~x+y
    unconSim = krigeSimCE1(z~1, newdata = loc, model = modVgm, n=tr)
  }else{
    coordinates(loc) = ~x+y
    g.dummy = gstat(formula = z~1, dummy = TRUE, beta = 0,
                     model = modVgm, nmax = 50)
    unconSim = predict(g.dummy, loc, nsim = tr)
  }
  return(unconSim)
}

predict.surfaceGrid<- function(object,x){
  interp.surface( object, x)
}
predict.multivariateSurfaceGrid<- function(object,x){
  dimZ<- dim( object$z)
  L<- dimZ[3]
  out<- matrix( NA, nrow= nrow(x), ncol=L)
  for ( l in 1:L){
    out[,l]<- interp.surface(
      list( x=object$x,y=object$y, z=object$z[,,l]) , x)
  }
  return( out)
}
predict.constantValue<- function(object, x){
  n<- length(object$values)
  m<- nrow( x)
  return( matrix( object$values, nrow=m, ncol=n, byrow=TRUE ) )
}

lk_generator = function(x, y, kappa, tausq){
  gridList =  list(x = x, y = y)
  xPoints =  make.surface.grid(gridList)
  rhoTemp =  8 + 2*pnorm(xPoints[,1] + 2*xPoints[,2], mean=.25, sd =.3)
  # image.plot(as.surface(xPoints, rhoTemp))
  rho.object =  as.surface(xPoints, rhoTemp)
  class(rho.object) =  "surfaceGrid"
  
  taper<- pnorm(xPoints[,1] - xPoints[,2],
                mean = 0, sd=.15)
  a.wghtTemp<- 5.1 + taper 
  a.wghtObjectA <- as.surface(xPoints, a.wghtTemp)
  class(a.wghtObjectA)<- "surfaceGrid"
  
  sDomain = apply(coords, 2, range)
  LKInfo = LKrigSetup(sDomain, NC = 10, NC.buffer = 0, nlevel = 3, 
                      a.wghtObject = a.wghtObjectA,
                      nu = kappa,
                      rho.object = rho.object,
                      sigma = tausq)
  
  z =  LKrig.sim(xPoints, LKInfo)
  return(z)
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
  
  z = st_generator(coords, pars, method = "CE", tr = 1)
  z = z[[1]]
  dat = cbind(as.matrix(coords), z)
  
  # Using fields pkg:
  # obj =  circulantEmbeddingSetup(grid = list(x=x, y=y),
  #                                cov.args=list( Covariance="Matern", aRange=phi, smoothness=kappa))
  # z = as.vector(sqrt(sigmasq)*circulantEmbedding(obj))
  # z = z + sqrt(tausq)*rnorm(n)
  
  z.g = data.frame(coords, data = z)
  coordinates(z.g) = ~x+y
  v.sample = gstat::variogram(data~1, z.g)
  v = mean(tail(v.sample)$gamma)
  r0 = max(v.sample$dist)/3
  # our vecchia
  loglik_true_mr = loglik_mr(pars, dat, m, K)
  # pars_mr = optim(par = c(v/2, r0, 1, v/2), fn = ng_loglik_mr,
  #                dat = dat, m = m, K = K,
  #                method = "L-BFGS-B",
  #                lower = c(v/10, r0/10, 0.1, v/10),
  #                upper = c(v*100, r0*10, 10, v*10 ))
  
  # GPvecchia pkg
  vecchia.approx = vecchia_specify(dat[,1:2], m = K, ordering = "maxmin", conditioning = "NN")
  loglik_true_gp = vecchia_likelihood(dat[,3], vecchia.approx, covparms = pars[1:3], nuggets = pars[4])
  
  # pars_gp = optim(par = c(v/2, r0, 1, v/2), fn = ng_loglik_gp,  
  #                 dat = dat, vecchia.approx = vecchia.approx, 
  #                 method = "L-BFGS-B", 
  #                 lower = c(v/10, r0/10, 0.1, v/10), 
  #                 upper = c(v*100, r0*10, 10, v*10 ))
  
  # out conservative method
  loglik_true_cons = loglik_cons(pars, dat, m, K) 
  # pars_cons = optim(par = c(v/2, r0, 1, v/2), fn = ng_loglik_cons,
  #                 dat = dat, m = m, K = K,
  #                 method = "L-BFGS-B",
  #                 lower = c(v/10, r0/10, 0.1, v/10),
  #                 upper = c(v*100, r0*10, 10, v*10 ))
  
  #true likelihood
  if(n < 10000){
    D = dist1(as.matrix(coords))
    V = MaternFun(D, c(sigmasq, phi, kappa)) + tausq*diag(rep(1, nrow(D)))
    loglik_true = loglik(z, V)
  }else{
    loglik_true = loglik_cons(pars, dat, m = 10, K = 40) 
  }
  
  # if(n < 5000){
  #   pars_mle = optim(par = c(v/2, r0, 1, v/2), fn = ng_loglik,
  #                   dat = dat, 
  #                   method = "L-BFGS-B",
  #                   lower = c(v/10, r0/10, 0.1, v/10),
  #                   upper = c(v*100, r0*10, 10, v*10 ))
  # }
  
  reletive_loglik = c(loglik_true_mr/loglik_true, 
                      loglik_true_cons/loglik_true,
                      loglik_true_gp/loglik_true)
  #bias_mr = pars_mr$par - pars
  #bias_gp = pars_gp$par - pars
  return(reletive_loglik)
}

num.core = detectCores()

sigmasq = 10
kappa = 2.25
tausq = 0.25
m = 5
K = 10
for(N in c(50, 100, 150)){
  res1 = list()
  res2 = list()
  res3 = list()
  for(i in 1:5){
    res1 = append(res1, sim(N, m, K, pars = c(sigmasq, 0.025, kappa, tausq)))
    
  }


  save.image(paste0("v", N, ".RData"))

}

  

