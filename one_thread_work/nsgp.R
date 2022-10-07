library(SpatialTools)
library(gstat)
library(sp)
library(MASS)

K1_sigma = function(coord, sigma){
  sigma[1]*(2+cos(sigma[2]*(coord[1]-sigma[3])))+cos(sigma[2]*(coord[2]-sigma[3]+0.1))
}
K1 = function(coord1, coord2, theta=0.5, lambda=c(0.025, 0.09), sigma=c(5,5,3), kappa=1.25, phi=1){
  #lambda is of length 2, sigma is of length 3
  coord1 = as.vector(as.matrix(coord1))
  coord2 = as.vector(as.matrix(coord2))
  rot = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),
                 nrow = 2, byrow = TRUE)
  G.inv = rot%*%diag(1/lambda)%*%t(rot)
  dist.G = t(coord1-coord2)%*%G.inv%*%(coord1-coord2)
  if(dist.G){
    m = 2^(1-kappa)/gamma(kappa)*(sqrt(2*kappa)*dist.G/phi)^kappa*
      besselK(sqrt(2*kappa)*dist.G/phi, kappa)
  }else{m = 1}
  
  K1_sigma(coord1, sigma)*K1_sigma(coord2, sigma)*m
}

K2 = function(coord1, coord2, sigmasq = 2, phi = 5, kappa = 1.25){
  coord1 = as.vector(as.matrix(coord1))
  coord2 = as.vector(as.matrix(coord2))
  S1 = diag(coord1^2)
  S2 = diag(coord2^2)
  S12 = (S1+S2)/2
  Q12 = t(coord1-coord2)%*%sqrt(solve(S12))%*%(coord1-coord2)
  r.Q12 = sqrt(Q12)
  if(r.Q12){
    m = 2^(1-kappa)/gamma(kappa)*(sqrt(2*kappa)*r.Q12/phi)^kappa*
    besselK(sqrt(2*kappa)*r.Q12/phi, kappa)
    }else{m = 1}
  
  sigmasq*det(S1)^(1/4)*det(S2)^(1/4)/det(S12)^(1/2)*m
}


K3 = function(coord1, coord2, sigmasq = 2, phi = 5, kappa = 1.25){
  coord1 = as.vector(as.matrix(coord1))
  coord2 = as.vector(as.matrix(coord2))
  d = length(coord1)
  l1 = abs(coord1[1])+abs(coord1[2])/2
  l2 = abs(coord2[1])+abs(coord2[2])/2
  d12 = sqrt(sum((coord1-coord2)^2))/sqrt((l1+l2)/2)
  if(d12){
    m = 2^(1-kappa)/gamma(kappa)*(sqrt(2*kappa)*d12/phi)^kappa*
      besselK(sqrt(2*kappa)*d12/phi, kappa)
  }else{m=1}
  
  2^(d/2)*l1^(d/4)*l2^(d/4)/((l1+l2)/2)^(d/2)*m
}

K4_pair = function(i, j, coords, l1, l2, g1, g2, sigmasq, phi, kappa){
  #given the processes and index i and j,
  #get the covariance between si and sj
  l1i = l1[i]
  l1j = l1[j]
  l2i = l2[i]
  l2j = l2[j]
  
  g1i = g1[i]
  g1j = g1[j]
  g2i = g2[i]
  g2j = g2[j]
  
  di = sqrt(l1i)
  dj = sqrt(l1j)
  
  Gi = matrix(c(g1i, -g2i, g2i, g1i), nrow = 2, byrow = TRUE)/di
  Li = diag(c(l1i, l2i))
  Si = Gi%*%Li%*%t(Gi)
  
  Gj = matrix(c(g1j, -g2j, g2j, g1j), nrow = 2, byrow = TRUE)/dj
  Lj = diag(c(l1j, l2j))
  Sj = Gj%*%Lj%*%t(Gj)
  
  Sij = (Si+Sj)/2
  eij = eigen(Sij)
  r.Sij = eij$vectors%*%diag(1/sqrt(eij$values))%*%t(eij$vectors)
    
  Q = t(coords[i,]-coords[j,])%*%r.Sij%*%(coords[i,]-coords[j,])
  r.Q = sqrt(Q)
  if(r.Q){
    m = 2^(1-kappa)/gamma(kappa)*(sqrt(2*kappa)*r.Q/phi)^kappa*
      besselK(sqrt(2*kappa)*r.Q/phi, kappa)
  }else{m=1}
  
  
  sigmasq*det(Si)^(1/4)*det(Sj)^(1/4)/det(Sij)^(1/2)*m
  
}
K4 = function(coords, sigmasq = 2, phi = 5, kappa = 1.25){
  #output covariance matrix 
  coords  = as.matrix(coords)
  n = nrow(coords)
  D = dist1(as.matrix(coords))
  #generate GP for log lambda2 process
  ll2 = rmvnorm(nsim = 1, mu = rep(0, n), V = 0.5*exp(-D/2))
  l2 = exp(ll2)
  #generate GP for gamma1 process
  g1 = rmvnorm(nsim = 1, mu = rep(0, n), V = 0.5*exp(-D/2))
  #generate GP for gamma2 process
  g2 = rmvnorm(nsim = 1, mu = rep(0, n), V = 0.5*exp(-D/2))
  #calculate lambda1
  l1 = g1^2+g2^2
  
  sapply(1:n,function(i){
    sapply(1:n,function(j){
      K4_pair(i, j, coords, l1, l2, g1, g2, sigmasq, phi, kappa)
      })
    }) 
  
}

N = 30
x = c(1:N)/(N/10)
y = c(1:N)/(N/10)
coords = expand.grid(x, y)
names(coords) = c("x", "y")

n = N*N

V1 = sapply(1:n,function(i){
  sapply(1:n,function(j){
    K1(coords[i,], coords[j,])
  })
})

V2 = sapply(1:n,function(i){
  sapply(1:n,function(j){
    K2(coords[i,], coords[j,])
  })
})

V3 = sapply(1:n,function(i){
  sapply(1:n,function(j){
    K3(coords[i,], coords[j,])
  })
})

V4 = K4(coords)


test.plot = function(n, coV, coords){
  z = rmvnorm(nsim = 1, mu = rep(0, n), V = coV)
  coordinates(coords) = ~x+y
  geo.z = SpatialPointsDataFrame(coords, as.data.frame(z))
  spplot(geo.z)
}

test.plot(n, V1, coords)
test.plot(n, V2, coords)
test.plot(n, V3, coords)
test.plot(n, V4, coords)



