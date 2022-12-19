
removeDiags <- function(X){
  out <- X
  for(j in 1:dim(X)[1]){
    diag(out[j,,]) <- 0
  }
  return(out)
}

Obj.Func <- function(Ts, Omega, lambda1, lambda2){
  out <- 0
  for(j in 1:dim(Ts)[1]){
    out <- out + sum((Ts[j,,] - matrix(diag(Omega[j,,]), nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=TRUE)  - matrix(diag(Omega[j,,]),nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=FALSE) + 2*Omega[j,,])^2)
  }
  return(out + lambda1*sum(abs(removeDiags(Omega))) + lambda2*sum(apply(removeDiags(Omega), c(2,3), function(x){sqrt(sum(x^2))})))
}


Loss.Func <- function(Ts, Omega){
  out <- 0
  for(j in 1:dim(Ts)[1]){
    out <- out + sum((Ts[j,,] - matrix(diag(Omega[j,,]), nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=TRUE)  - matrix(diag(Omega[j,,]),nrow=dim(Omega)[2], ncol =dim(Omega)[2], byrow=FALSE) + 2*Omega[j,,])^2)
  }
  return(out)
}



prox.function <- function(Y, lambda1, lambda2){
  temp <- pmax(abs(Y) - lambda1, 0)*sign(Y)
  out <- array(0, dim = dim(Y))
  for(kk in 1:(dim(Y)[2]-1)){
    for(jj in (kk+1):dim(Y)[3]){
      if(kk != jj){
        t0 <- sqrt(sum(temp[,kk,jj]^2))
        out[,jj,kk] <- max(1 - lambda2/t0, 0)*temp[,kk,jj]
      } 
    }
  }
  for(j in 1:dim(Y)[1]){
    out[j,,] <- out[j,,] + t(out[j,,])
    diag(out[j,,]) <- diag(Y[j,,])
  }
  return(out)
}




grad.loss <- function(Ts, Omega){
  out <- array(0, dim=dim(Ts))
  for(j in 1:dim(Ts)[1]){
    for(k in 1:dim(Ts)[2]){
      for(l in k:dim(Ts)[3]){
        if(k == l){
          for(mm in c(1:dim(Ts)[2])[-k]){
            out[j,k,l] <- out[j,k,l] + 4*Omega[j,k,k] - 4*(Ts[j,k,mm] - Omega[j,mm,mm] + 2*Omega[j,k,mm])
          }
        } else {
          out[j,k,l] <- 4*(Ts[j,k,l] - Omega[j,k,k] - Omega[j,l,l]) + 8*Omega[j,k,l]
        }  
      }
    }
    out[j,,] <- out[j,,] + t(out[j,,])
    diag(out[j,,]) <- diag(out[j,,])/2
  }
  return(out)
}


# ------------------------------
# function inputs 
# ------------------------------
# alpha0 <- 1
# tau <- 1e-4
# max.iter <- 1000
# lambda1 <- .9
# lambda2 <- 1
# tol <- 1e-8
# X.init <- Omega
# U.init <- array(0, dim=dim(Omega))
# Z.init <- Omega
# -------------------------------

MCIScompute <- function(Ts,
  X.init = NULL,
  Z.init = NULL,
  U.init = NULL,
  alpha0 = 1, 
  tau = 1e-4,
  max.iter = 5e3,
  lambda1,
  lambda2,
  tol.obj = 1e-8,
  tol.diff = 1e-8,
  very.quiet = TRUE
) {

  alpha <- alpha0
  X.new <- X.init
  Z.new <- Z.init
  U.new <- U.init
  loss.prev <- Loss.Func(Ts, Z.new)
  obj.func <- rep(0, max.iter)

  for(iteration in 1:1000){
    
    grad.temp <- grad.loss(Ts, Z.new)
    lineSearch <- TRUE

    while(lineSearch){
      X.temp <- Z.new - alpha*U.new - alpha*grad.temp
      X.new <- prox.function(X.temp, alpha*lambda1, alpha*lambda2)
      loss.temp <- Loss.Func(Ts, X.new)
      if(loss.temp <= loss.prev + sum(grad.temp*(X.new - Z.new)) + sum((X.new - Z.new)^2)/(2*alpha)){
        lineSearch <- FALSE
      } else {
        alpha <- alpha*.5
      }
    }
    
    
    Z.new <- X.new + alpha*U.new
    for(kk in 1:dim(Z.new)[1]){
      eo <- eigen(Z.new[kk,,])
      if(!all(eo$val > tau)){
        Z.new[kk,,] <- crossprod(t(eo$vec), pmax(eo$val, tau)*t(eo$vec))
        #no.Need <- FALSE
      } 
    }

    U.temp <- array(0, dim=dim(Z.new))
    for(kk in 1:dim(Z.new)[1]){
      U.temp[kk,,] <- U.new[kk,,] + (X.new[kk,,] - Z.new[kk,,])/alpha
      U.temp[kk,,] <- (U.temp[kk,,] + t(U.temp[kk,,]))/2
    }

    U.new <- U.temp
    loss.prev <- Loss.Func(Ts, Z.new)
    #alpha <- min(1.1*alpha, alpha0)

    obj.func[iteration] <- loss.prev + lambda1*sum(abs(removeDiags(Z.new))) + lambda2*sum(apply(removeDiags(Z.new), c(2,3), function(x){sqrt(sum(x^2))}))
    if(!very.quiet){
      cat(obj.func[iteration], "\n")
    }
    if(iteration > 5){
      if(all(abs(obj.func[iteration:(iteration-2)] - obj.func[(iteration-1):(iteration-3)]) < tol.obj*abs(obj.func[iteration]))){
        if(max(abs(X.new - Z.new)) < tol.diff){
          break
        }
      }
    }
  }

  return(list("U.new" = U.new, "Z.new" = Z.new, "X.new" = X.new))
}






# out <- MCIScompute(Ts,
#   X.init = Omega,
#   Z.init = Omega,
#   U.init = array(0, dim= c(J, p, p)),
#   alpha0 = 1, 
#   tau = 1e-4,
#   max.iter = 5e3,
#   lambda1 = 3,
#   lambda2 = 1,
#   tol.obj = 1e-10,
#   tol.diff = 1e-8,
#   very.quiet = FALSE
# )





# T <- Ts[1,,]

# ----------------------
# diag object funtion is 
## sum sum_{i \neq j} (T_{ij} - w_{ii} - w_jj)^2
# ----------------------

MCISpath <- function(Ts, TsVal, 
      nlambda1 = 25, nlambda2 = 25, 
      deltalambda1 = 0.01, 
      deltalambda2 = 0.01, 
      alpha0 = 1,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-6,
      tol.diff = 1e-6,
      quiet = FALSE,
      very.quiet = TRUE){


  getDiagT <- function(T, max.iter = 1e4){
    p <- dim(T)[1]
    omega <- rep(1, dim(T)[1])
    omega.prev <- omega
    for(kk in 1:max.iter){
      for(jj in 1:p){
        omega[jj] <- sum(T[jj,-jj] - omega[-jj])/(p-1)
      }
      if(sum(abs(omega - omega.prev)) < 1e-12){
        break
      } else {
        omega.prev <- omega
      }
    }
    return(omega)
  }

  getDiagTs <- function(Ts){
    Omega <- array(0, dim(Ts))
    for(kk in 1:dim(Ts)[1]){
      diag(Omega[kk,,]) <- getDiagT(Ts[kk,,])
    }
    return(Omega)
  }


  # ---------------------------------------
  # get candidate tuning parameter pairs 
  # ---------------------------------------
  # nlambda1 <- 25
  # deltalambda1 <- .01
  # nlambda2 <- 25
  # deltalambda2 <- .01

  out <- getDiagTs(Ts)
  t0 <- grad.loss(Ts, out)
  lambda2max <- max(apply(t0, c(2,3), function(x){sqrt(sum(x^2))}))
  lambda2seq <- 10^seq(log10(lambda2max), log10(deltalambda2*lambda2max), length=nlambda2+2)[-c(1,2)]
  lambda1mat <- matrix(0, nrow=nlambda2, ncol=nlambda1)
  lambda1max.candidate <- max(removeDiags(abs(t0)))
  lambda1.candidates <- 10^seq(log10(lambda1max.candidate), log10(.01*lambda1max.candidate), length=100)
  U.init <- array(0, dim= dim(Ts))
  X.init <- Ts
  Z.init <- Ts

  for(kk in 1:nlambda2){

    for(jj in 1:100){
       out <- MCIScompute(Ts,
        X.init = X.init,
        Z.init = Z.init,
        U.init = U.init,
        alpha0 = alpha0, 
        tau = tau,
        max.iter = max.iter,
        lambda1 = lambda1.candidates[jj],
        lambda2 = lambda2seq[kk],
        tol.obj = tol.obj,
        tol.diff = tol.diff,
        very.quiet = very.quiet
      )

      X.init <- out$X.new
      Z.init <- out$Z.new
      U.init <- out$U.new
      if(sum(removeDiags(X.init)!=0) > 0){
        break
      } else {
        cat(jj, "\n")
      }
    }
    if(jj > 1){
      lambda1mat[kk,] <- 10^seq(log10(lambda1.candidates[jj-1]), log10(deltalambda1*lambda1.candidates[jj-1]), length=nlambda1)
    } else {
      lambda1mat[kk,] <- 10^seq(log10(lambda1.candidates[jj]), log10(deltalambda1*lambda1.candidates[jj]), length=nlambda1)
    }
  }

  Omega.mat <- array(0, dim=c(dim(Ts)[1], dim(Ts)[2], dim(Ts)[3], nlambda1, nlambda2))
  results <- matrix(Inf, nrow=nlambda1, ncol=nlambda2)

  for(kk in 1:nlambda2){
    for(jj in 1:nlambda1){
      out <- MCIScompute(Ts,
        X.init = X.init,
        Z.init = Z.init,
        U.init = U.init,
        alpha0 = alpha0, 
        tau = tau,
        max.iter = max.iter,
        lambda1 = lambda1mat[kk,jj],
        lambda2 = lambda2seq[kk],
        tol.obj = tol.obj,
        tol.diff = tol.diff,
        very.quiet = very.quiet)

      X.init <- out$X.new
      Z.init <- out$Z.new
      U.init <- out$U.new

      Omega.mat[,,,jj,kk] <- out$X.new
      results[jj,kk] <- Loss.Func(TsVal, Omega.mat[,,,jj,kk])
      if(jj > 4){
        if(all(results[jj:(jj-3),kk] > results[(jj-1):(jj-4),kk])){
          break
        }
      }
      if(!quiet){
        cat(results[jj,kk], "\n")
      }
    }
  }

  return(list("Omega" = Omega.mat, 
              "results" = results,
              "lambda1" = lambda1mat,
              "lambda2" = lambda2seq))
}


# t0 <- MCISpath(Ts = Ts, TsVal = TsVal, 
#       nlambda1 = 10, 
#       nlambda2 = 10, 
#       deltalambda1 = 0.01, 
#       deltalambda2 = 0.01, 
#       alpha0 = 1,  
#       max.iter = 5e3,
#       tol.obj = 1e-6,
#       tol.diff = 1e-6,
#       tau = 1e-4,
#       quiet = FALSE,
#       very.quiet = TRUE)
