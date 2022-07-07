


remove(list = ls())
library(maxLik)
library(survival)
library(foreach)
library(doParallel)
library(MASS)


seed = 123
N = 10000  # total sample size
n = 200  #  the size of the i dataset
B = N / n # the Number of sub datasets,  N / n
p = 2
repetition = 500


beta.true = c(-1,1)
gamma.true = 0.5
theta.true = exp(gamma.true)
p = length(beta.true)
a0 = 0.3  # adjust cure rate, low with low, low cure = 0.3, high cure = 1.25
c0 = 0.4 # adjust censoring rate, larger c0, then higher censoring rate

F0 = function(t){
  f = t*ifelse(0<=t&t<=1,1,0) + 1*ifelse(t>1,1,0)
  return(as.vector(f))
}
g = function(beta,x){
  return(as.vector(exp(x%*%beta) )) # dim is n*1
}
dg = function(beta,x){
  return(x * g(beta,x)) # x*g(beta,x), dim is n*p, x is n*p
}
ddg = function(beta,x){
  return(sweep(apply(x, 1, function(u) u%*%t(u)), MARGIN = 2, g(beta,x), '*') ) # xx * g(beta.x), dim is (p^2)*n, the j col is as.vector(xx) of j person
}
Q = function(g.value, Delta, W){
  n = length(g.value)
  return( (cumsum((W*g.value*Delta)[n:1])[n:1] + sum(W*g.value*(1-Delta)))/n )  # Q(Y.i), dim is n*1
}
Likhd = function(beta, data, W){
  ## log Likelihood, dim is 1*1, f = sum W*delta*log(g/Q)
  x = as.matrix(data[,-(1:2)])  # dim is n*p
  tau = max(data$Y[data$delta==1])
  Delta = as.numeric(data$Y <= tau)  # Delta=I(Y<=tau)
  g.value = g(beta,x)  # dim is n*1, g(beta,x)
  Q.value = Q(g.value,Delta,W)  # dim is n*1, Q(Y.i)
  g.value[(data$delta==0)] = 1
  Q.value[(data$delta==0)] = 1
  f = sum(W*data$delta*log(g.value/Q.value), na.rm = TRUE)
  return(f)
}
dLikhd = function(beta, data, W){
  # socre function, dim is p*1, f = sum W*delta*[g'/g - Q'/Q]
  x = as.matrix(data[,-(1:2)])  # dim is n*p
  tau = max(data$Y[data$delta==1])
  Delta = as.numeric(data$Y <= tau)  # Delta=I(Y<=tau)
  g.value = g(beta,x)  # dim is n*1, g(beta,x)
  Q.value = Q(g.value,Delta,W)  # dim is n*1, Q(Y.i)
  dg.value = dg(beta,x)  # dim is n*p, x*g(beta,x)
  dQ.value = apply(dg.value, 2, Q, Delta, W)  # dim is n*p  
  g.value[data$delta==0] = 1
  Q.value[data$delta==0] = 1
  f = apply(W*data$delta*(dg.value/g.value - dQ.value/Q.value), 2, sum, na.rm = TRUE)  # dim is p*1
  return(as.vector(f))
}
ddLikhd = function(beta, data, W){
  # second derivative of log-Likelihood, dim is p*p, f = sum W*delta*[(g''g-g'g')/gg - (Q''Q-Q'Q')/QQ]
  p = length(beta)
  x = as.matrix(data[,-(1:2)])  # dim is n*p
  tau = max(data$Y[data$delta==1])
  Delta = as.numeric(data$Y <= tau)  # Delta=I(Y<=tau)
  
  g.value = g(beta,x)  # dim is n*1, g(beta,x)
  dg.value = dg(beta,x)  # dim is n*p, g' = x*g(beta,x)
  ddg.value = ddg(beta,x)  # g''= xx*g(beta,x), dim is (p^2)*n
  
  Q.value = Q(g.value,Delta,W)  # dim is n*1, Q(Y.i)
  dQ.value = apply(dg.value, 2, Q, Delta, W)  # dim is n*p, derivative of Q(Y.i)
  ddQ.value = t(apply(ddg.value, 1, Q, Delta, W))  # dim is (p^2)*n, second derivative of derivative of Q(Y.i)
  
  ddg_g = sweep(ddg.value, 2, g.value, '*')  # g''(beta,x) * g(beta,x), dim is (p^2)*n
  dg_dg = apply(dg.value, 1, function(u) u%*%t(u))  # g'(beta,x)*g'(beta,x), dim is (p^2)*n
  g.value[(data$delta)==0] = 1
  Q.value[(data$delta==0)] = 1
  a1 = sweep(ddg_g-dg_dg, 2, g.value*g.value, '/')  # (g''g-g'g')/gg, dim is (p^2)*n
  
  ddQ_Q = sweep(ddQ.value, 2, Q.value, '*') # Q''*Q, dim is (p^2)*n
  dQ_dQ = apply(dQ.value, 1, function(u) u%*%t(u))  # Q'(beta,x)*Q'(beta,x), dim is (p^2)*n
  a2 = sweep(ddQ_Q-dQ_dQ, 2, Q.value*Q.value, '/')  # (Q''Q-Q'Q')/QQ, dim is (p^2)*n
  
  f = sweep(a1-a2, 2, W*data$delta,'*')  # W*delta*[(g''g-g'g')/gg - (Q''Q-Q'Q')/QQ], dim is (p^2)*n
  f = apply(f, 1, sum)  # dim is (p^2)*1
  f = matrix(f, nrow = p, ncol = p, byrow = FALSE)
  return(f)
}
ddLikhd2 = function(beta, data, W){
  # second derivative of log-Likelihood, dim is p*p, f = sum W*delta*[(g''g-g'g')/gg - (Q''Q-Q'Q')/QQ]
  p = length(beta)
  x = as.matrix(data[,-(1:2)])  # dim is n*p
  tau = max(data$Y[data$delta==1])
  Delta = as.numeric(data$Y <= tau)  # Delta=I(Y<=tau)
  g.value = g(beta,x)  # dim is n*1, g(beta,x)
  g.value.record = g.value  # dim is n*1, g(beta,x)
  dg.value = dg(beta,x)  # dim is n*p, g' = x*g(beta,x)
  Q.value = Q(g.value,Delta,W)  # dim is n*1, Q(Y.i)
  dQ.value = apply(dg.value, 2, Q, Delta, W)  # dim is n*p, derivative of Q(Y.i)
  g.value[(data$delta)==0] = 1
  Q.value[(data$delta==0)] = 1
  f = matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:i) {
      ddg.ij = x[,i]*x[,j]*g.value.record  # dim is n*1
      ddQ.ij = Q(ddg.ij, Delta, W)  # dim is n*1, second derivative of derivative of Q(Y.i)
      f.ij = W * data$delta * ((ddg.ij*g.value - dg.value[,i]*dg.value[,j])/g.value^2 - (ddQ.ij*Q.value - dQ.value[,i]*dQ.value[,j])/Q.value^2)  # W*delta*[(g''g-g'g')/gg - (Q''Q-Q'Q')/QQ], dim is n*1  
      f[i,j] = sum(f.ij, na.rm = TRUE)  # dim is 1*1
    }
  }
  f[upper.tri(f)] = f[lower.tri(f)]
  return(f)
}
J = function(beta, data, W){
  # Negative second derivative of log-Likelihood, dim is p*p, f = sum W*delta*[(g''g-g'g')/gg - (Q''Q-Q'Q')/QQ]
  p = length(beta)
  n = nrow(data)
  x = as.matrix(data[,-(1:2)])  # dim is n*p
  tau = max(data$Y[data$delta==1])
  Delta = as.numeric(data$Y <= tau)  # Delta=I(Y<=tau)
  
  g.value = g(beta,x)  # dim is n*1, g(beta,x)
  dg.value = dg(beta,x)  # dim is n*p, g' = x*g(beta,x)
  ddg.value = ddg(beta,x)  # g''= xx*g(beta,x), dim is (p^2)*n
  
  Q.value = Q(g.value,Delta,W)  # dim is n*1, Q(Y.i)
  dQ.value = apply(dg.value, 2, Q, Delta, W)  # dim is n*p, derivative of Q(Y.i)
  ddQ.value = t(apply(ddg.value, 1, Q, Delta, W))  # dim is (p^2)*n, second derivative of derivative of Q(Y.i)
  
  ddg_g = sweep(ddg.value, 2, g.value, '*')  # g''(beta,x) * g(beta,x), dim is (p^2)*n
  dg_dg = apply(dg.value, 1, function(u) u%*%t(u))  # g'(beta,x)*g'(beta,x), dim is (p^2)*n
  g.value[(data$delta)==0] = 1
  Q.value[(data$delta==0)] = 1
  a1 = sweep(ddg_g-dg_dg, 2, g.value*g.value, '/')  # (g''g-g'g')/gg, dim is (p^2)*n
  
  ddQ_Q = sweep(ddQ.value, 2, Q.value, '*') # Q''*Q, dim is (p^2)*n
  dQ_dQ = apply(dQ.value, 1, function(u) u%*%t(u))  # Q'(beta,x)*Q'(beta,x), dim is (p^2)*n
  a2 = sweep(ddQ_Q-dQ_dQ, 2, Q.value*Q.value, '/')  # (Q''Q-Q'Q')/QQ, dim is (p^2)*n
  
  f = sweep(a1-a2, 2, W*data$delta,'*')  # W*delta*[(g''g-g'g')/gg - (Q''Q-Q'Q')/QQ], dim is (p^2)*n
  f = apply(f, 1, sum)  # dim is (p^2)*1
  f = matrix(f, nrow = p, ncol = p, byrow = FALSE)
  f = -f
  return(f)
}
Transform = function(x) {
  r.num = nrow(x) 
  y = list()
  for(r in 1:r.num){
    y[[r]] = x[r,]
  }
  return(y)
}
Get.data = function(n){
  # n: sample size
  
  x1 = runif(n,a0,a0+1)
  # x2 = rnorm(n, mean = a0, sd = 1/12)
  x2 = rbinom(n,1,0.5)
  x = cbind(x1,x2)
  
  cure = exp(-g(beta.true,x)*theta.true) # cure probabilities
  is.cure = rbinom(n, 1, cure)
  T0 = NULL 
  u <- runif(n, min = 0, max = 1-cure)
  for (i in 1:n){
    if(is.cure[i]==0){
      T0[i] = uniroot(function(t) {1 - exp(-g(beta.true,x)[i]*theta.true*F0(t)) - u[i]}, c(0, 100))$root  # uncured survival time
    }else{
      T0[i] = 10000  # cured survival time
    }
  }
  
  C <- rexp(n, c0)  # censoring time
  Y = pmin(T0, C)  # observed time
  delta = as.numeric(T0<=C)  # censoring indicator
  data = data.frame(Y,delta,x1,x2)
  
  cureR = mean(cure)  # cure rate
  cenR = sum(delta==0)/n  # censoring rate
  
  return(list(data=data, cenR=cenR, cureR=cureR))
}
Quadraticform = function(b,A){
  # get b'A^-1b, where A is a positive definite matrix
  L = t(chol(A))  #  Choleski Decomposition, get a lower triangular matri
  x = backsolve(L, b, upper.tri = FALSE, transpose = FALSE) # Solve an Lower Triangular System
  f = sum(x^2)
  return(f)
}
Getsummary1 = function(){
  # summary for beta
  beta.est = beta.se = matrix(0, nrow = repetition, ncol = p)
  cureR = cenR =  rep(0, repetition)
  for(i in 1:repetition){
    beta.est[i,] = MyRsult[[i]]$beta.est
    beta.se[i,] = MyRsult[[i]]$beta.se
    cureR[i] = MyRsult[[i]]$cureR; cenR[i] = MyRsult[[i]]$cenR
  }
  ccureR = mean(cureR)  # cure rate
  ccenR = mean(cenR)  # censoring rate
  
  beta.est.avrge = apply(beta.est, 2, mean)
  beta.SD = apply(beta.est, 2, sd)
  beta.SE = apply(beta.se, 2, mean)
  beta.bias = beta.est.avrge - beta.true
  beta.Abias = apply(t(apply(beta.est, 1, function(x) abs(x-beta.true))), 2, mean)
  
  get.CP = function(TrueVaule, est, se){
    # TrueVaule: dim is p*1
    # est: dim is repetition*P
    # se: dim is repetition*p
    p = length(TrueVaule)
    if(!is.matrix(se)){ se = matrix(se, nrow = repetition, ncol = p, byrow = TRUE)}  # if input SD, make SD to be a matrix
    TrueVaule.mat = matrix(TrueVaule, nrow = repetition, ncol = p, byrow = TRUE)
    CP = apply((est - qnorm(0.975,0,1)*se)<=TrueVaule.mat & TrueVaule.mat<=(est + qnorm(0.975,0,1)*se), 2, mean)
    return(CP)
  }
  CP.beta = get.CP(beta.true, beta.est, beta.se)
  
  # summary 
  runningtime = difftime(Sys.time(), timestart, units="mins")
  runningtime = round(runningtime, 2)
  our = cbind(beta.bias, beta.SD, beta.SE, CP.beta)
  our = cbind(our, runningtime)
  rownames(our) = c('beta1.Renew', 'beta2.Renew')
  colnames(our) = c('Bias', 'SD', 'SE', 'CP', 'Time')
  
  
  cat('===========================================================', '\n',
      'Run Time                :', runningtime, 'minus', '\n',
      'Nample size    n        :', n,   '\n',
      'Number batach  B        :', B,   '\n',
      'Right-censoring rate    :', ccenR, '\n',
      'Cure rate               :', ccureR, '\n'
  )
  print(our)
  cat('===========================================================', '\n','\n','\n')
  return(our)
}
Renew.est1 = function(r){
  # get beta.est and beta.se  
  set.seed(seed + r)
  tol = 1e-5;
  max_iter = 100;
  Jsum = diag(0,p,p)  # J_1(beta.1)+...+J_b(beta.b)
  betahat = rep(0,p)  # initial value
  cureR = cenR = NULL
  for (b in 1:B){
    dat = Get.data(n)
    data = dat$data
    cureR[b] = dat$cureR
    cenR[b] = dat$cenR
    n = nrow(data); p = ncol(data)-2
    data = data[order(data[,1]), ]  # increasing order by Y
    x = as.matrix(data[,-(1:2)])
    tau = max(data$Y[data$delta==1])
    Delta = as.numeric(data$Y <= tau)  # Delta=I(Y<=tau)
    W = rep(1,n)  # disturbed value
    
    # get beta.est
    betahat_old = betahat # beta.b-1, also is the initial value beta.b^1
    J.b = J(beta = betahat_old, data = data, W=W)  # J.b(beta.b-1)
    JJ = chol(Jsum + J.b)  # Choleski factorization of Jsum.b-1 + J.b(beta.b-1), denoted as JJ
    L = t(JJ)
    
    for (r in 1:max_iter){
      U.b = dLikhd(beta = betahat, data = data, W=W)  # U.b(beta.b^r) 
      JU = drop(Jsum %*% (betahat_old-betahat)) # Jsum.b-1 %*% (beta.b-1 - beta.b^r)
      UU = JU + U.b  # inceremental equation: Jsum.b-1 %*% (beta.b-1 - beta.b^r) + U.b(beta.b^r), denoted as UU
      d_beta=backsolve(JJ,forwardsolve(L,UU))  # JJ^-1 %*% UU, that is using  Choleski factorization to solve X: JJ*X=UU
      eps = sqrt(sum(UU^2))
      if(eps<tol){
        break
      }else {
        betahat=betahat+d_beta;
      }
    }
    Jnew = J(beta = betahat, data = data, W = W)  # J.b(beta.b)
    Jsum = Jsum + Jnew  # J_1(beta.1)+...+J_b(beta.b)
  }
  
  beta.est = betahat
  beta.se = sqrt(diag(solve(Jsum)))  # get beta.se
  cureR = mean(cureR)
  cenR = mean(cenR)
  
  return(list(beta.est=beta.est, beta.se=beta.se
              ,cureR=cureR, cenR=cenR))
}


pklist = c('MASS','survival','maxLik')
varlist = c('p','n')

# get beta.est and beta.se
timestart<-Sys.time()
registerDoParallel(makeCluster(detectCores(logical = FALSE)))  # parallel
MyRsult = foreach(r = 1:repetition, .combine = rbind, .multicombine = TRUE, .packages = pklist, .export=varlist, .inorder=FALSE)  %dopar%  Renew.est1(r)
stopImplicitCluster()  
MyRsult = Transform(MyRsult)
MyRsult.beta = Getsummary1()







