# -----------------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# ------------------------------------------------------
source("~/blue/CompCovariance/Simulations/Functions/MCISpath.R")
source("~/blue/CompCovariance/Simulations/Functions/MCIS0path.R")
source("~/blue/CompCovariance/Simulations/Functions/COATpath.R")
source("~/blue/CompCovariance/Simulations/Functions/cclassoCV.R")
# ------------------------------------------------------

set.seed(uu)
nreps <- 50
params <- expand.grid(n = c(50, 100, 150), 
	p = rep(c(40, 80, 120, 160, 200), each=nreps), 
	model = c(1, 2, 3))
n <- params[uu,1]
p <- params[uu,2]
Model <- params[uu,3]
savename <- paste("~/blue/CompCovariance/Simulations/Results/Model", 
	Model, "/n", n, "_p",p, "_uu", uu%%nreps + 1 ,".RDS", sep="")

J <- 4
nval <- n

if(Model == 1){ 

	SigmaArray <- array(0, dim=c(J,p,p))
	for(ll in 1:2){
		for(jj in 1:p){
			for(kk in 1:p){
				if(abs(jj - kk) < 3){
					SigmaArray[ll,jj,kk] <- 0.3
				}
			}
		}
		diag(SigmaArray[ll,,]) <- 1
	} 
	for(ll in 3:4){
		for(jj in 1:p){
			for(kk in 1:p){
				if(abs(jj - kk) < 3){
					SigmaArray[ll,jj,kk] <- -0.2
				}
			}
		}
		diag(SigmaArray[ll,,]) <- 1
	} 

}

if(Model == 2){ 

	SigmaArray <- array(0, dim=c(J,p,p))
	for(jj in 1:p){
		for(kk in 1:p){
			if(abs(jj - kk) < p/4){
				SigmaArray[,jj,kk] <- 0.8^abs(jj-kk)
			}
		}
	}

	a <- 1
	b <- (p/4)
	SigmaArray[1,a:b,] <- SigmaArray[1,,a:b] <- 0
	SigmaArray[1,a:b, a:b] <- diag(1, p/4)

	a <- (p/4)+1
	b <- 2*(p/4)
	SigmaArray[2,a:b,] <- SigmaArray[2,,a:b] <- 0
	SigmaArray[2,a:b, a:b] <- diag(1, p/4)

	a <- 2*(p/4)+1
	b <- 3*(p/4)
	SigmaArray[3,a:b,] <- SigmaArray[3,,a:b] <- 0
	SigmaArray[3,a:b, a:b] <- diag(1, p/4)

	a <- 3*(p/4)+1
	b <- 4*(p/4)
	SigmaArray[4,a:b,] <- SigmaArray[4,,a:b] <- 0
	SigmaArray[4,a:b, a:b] <- diag(1, p/4)

}

if(Model == 3){

	SigmaArray <- array(0, dim=c(J,p,p))

	a <- 1
	b <- (p/2)
	for(kk in a:b){
		for(jj in a:b){
			SigmaArray[1,kk,jj] <- 0.9^(abs(kk-jj))
		}
	}

	a <- (p/6) + 1
	b <- (p/2) + (p/6)
	for(kk in a:b){
		for(jj in a:b){
			SigmaArray[2,kk,jj] <- 0.9^(abs(kk-jj))
		}
	}

	a <- 2*(p/6) + 1
	b <- (p/2) + 2*(p/6)
	for(kk in a:b){
		for(jj in a:b){
			SigmaArray[3,kk,jj] <- 0.9^(abs(kk-jj))
		}
	}

	a <- 3*(p/6) + 1
	b <- (p/2) + 3*(p/6)
	for(kk in a:b){
		for(jj in a:b){
			SigmaArray[4,kk,jj] <- 0.9^(abs(kk-jj))
		}
	}

	diag(SigmaArray[1,,]) <- diag(SigmaArray[2,,]) <- diag(SigmaArray[3,,]) <-  diag(SigmaArray[4,,]) <- 1
	for(kk in 1:4){
		SigmaArray[kk,,] <- diag(sqrt(seq(3, 1, length=p)))%*%SigmaArray[kk,,]%*%diag(sqrt(seq(3, 1, length=p)))
	}
}



compResults <- function(inputCov, SigmaArray, method){
	J <- dim(inputCov)[1]
	frobErrCov <- rep(0, J)
	specErrCov <- rep(0, J)
	frobErrCor <- rep(0, J)
	specErrCor <- rep(0, J)
	L1ErrCov <- rep(0, J)
	L1ErrCor <- rep(0, J)
	TPR <- rep(0, J)
	TNR <- rep(0, J)
	for(j in 1:J){
		frobErrCov[j] <- sqrt(sum((inputCov[j,,] - SigmaArray[j,,])^2))
		specErrCov[j] <- max(svd(inputCov[j,,] - SigmaArray[j,,])$d)
		L1ErrCov[j] <- max(apply(inputCov[j,,] - SigmaArray[j,,], 1, function(x){sum(abs(x))}))
		inputCor <- cov2cor(inputCov[j,,])
		trueCor <- cov2cor(SigmaArray[j,,])
		frobErrCor[j] <- sqrt(sum((inputCor - trueCor)^2))
		specErrCor[j] <- max(svd(inputCor - trueCor)$d)
		L1ErrCor[j] <- max(apply(inputCor - trueCor, 1, function(x){sum(abs(x))}))
		TPR[j] <- sum(inputCor!=0 & trueCor!=0)/sum(trueCor!=0)
		TNR[j] <- sum(inputCor==0 & trueCor==0)/sum(trueCor==0)
	}
	result.list <- list("frobErrCov" = frobErrCov,
		"specErrCov" = specErrCov,
		"frobErrCor" = frobErrCor, 
		"specErrCor" = specErrCor,
		"L1ErrCov" = L1ErrCov,
		"L1ErrCor" = L1ErrCor,
		"TPR" = TPR,
		"TNR" = TNR)
	
	names(result.list) <- paste(method, ".", names(result.list), sep="")
	return(result.list)
}

# -------------------------------------
# Generate data
# -------------------------------------
X <- array(0, dim=c(J, n, p))
Xval <- array(0, dim=c(J, n, p))
Y <- array(0, dim=c(J, n, p))
Yval <- array(0, dim=c(J, n, p))

for(j in 1:J){
  
	eo <- eigen(SigmaArray[j,,])
	SigmaXsqrt <- eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
	Y[j,,] <- matrix(rnorm(n*p), nrow=n)%*%SigmaXsqrt
	W <- exp(Y[j,,])
	X[j,,] <- W/rowSums(W)

	Yval[j,,] <- matrix(rnorm(n*p), nrow=n)%*%SigmaXsqrt
	W <- exp(Yval[j,,])
	Xval[j,,] <- W/rowSums(W)
}


# ---------------------------------------
# COAT
# ---------------------------------------
coatSoft <- list()
coatHard <- list()
coatSoftEst <- array(0, dim=c(J, p, p))
coatHardEst <- array(0, dim=c(J, p, p))

for(j in 1:J){
	coatSoft[[j]] <- coatPath(X[j,,],  Xval[j,,], soft = 1)
	coatSoftEst[j,,] <- coatSoft[[j]]$sigma
	coatHard[[j]] <- coatPath(X[j,,],  Xval[j,,], soft = 0)
	coatHardEst[j,,] <- coatHard[[j]]$sigma
}

Results <- compResults(coatSoftEst, SigmaArray, "coatSoft")
Results <- append(Results, compResults(coatHardEst, SigmaArray, "coatHard"))


# --------------------------------------
# MCIS
# --------------------------------------
Ts <- array(0, dim=c(J, p, p))
Tsval <- array(0, dim=c(J, p, p))
for(kk in 1:J){
	gMat <- matrix(0, p, p)
	for(jj in 1:p){
	    for(ll in 1:p){
	       gMat[jj,ll] <- sum(log(X[kk,,jj]/X[kk,,ll]))/n
	    }
	}
	for(jj in 1:p){
	    for(ll in 1:p){
	      Ts[kk,jj,ll] <- sum((log(X[kk,,jj]/X[kk,,ll]) - gMat[jj,ll])^2)/n
	    }
	}

	gMat <- matrix(0, p, p)
	for(jj in 1:p){
	    for(ll in 1:p){
	       gMat[jj,ll] <- sum(log(Xval[kk,,jj]/Xval[kk,,ll]))/n
	    }
	}
	for(jj in 1:p){
	    for(ll in 1:p){
	      Tsval[kk,jj,ll] <- sum((log(Xval[kk,,jj]/Xval[kk,,ll]) - gMat[jj,ll])^2)/n
	    }
	}


}


t0 <- MCISpath(Ts, Tsval, 
      nlambda1 = 20, nlambda2 = 10, 
      deltalambda1 = 0.0001, 
      deltalambda2 = 0.05, 
      alpha0 = 1,  
      tau = 1e-4,
      max.iter = 5e3,
      tol.obj = 1e-9,
      tol.diff = 1e-9,
      quiet = FALSE,
      very.quiet = TRUE)


getInd <- which(t0$results == min(t0$results), arr.ind=TRUE)
MCIS <- t0$Omega[,,,getInd[1,1],getInd[1,2]]
Results <- append(Results, compResults(MCIS, SigmaArray, "MCIS"))
ValErrs <- t0$results
Results <- append(Results, list("ValErrs" = ValErrs))

# ------------------------------------------
#
# ------------------------------------------
MCIS0 <- array(0, dim=c(J,p,p))
for(kk in 1:J){
  t0 <- MCIS0path(Ts[kk,,,drop=FALSE], Tsval[kk,,,drop=FALSE], 
                 nlambda1 = 25, 
                 deltalambda1 = 0.05, 
                 alpha0 = 1,  
                 tau = 1e-4,
                 max.iter = 5e3,
                 tol.obj = 1e-9,
                 tol.diff = 1e-9,
                 quiet = FALSE,
                 very.quiet = TRUE)
  getInd <- which(t0$results == min(t0$results))
  MCIS0[kk,,] <- t0$Omega[,,,getInd]
  ValErrs0 <- t0$results
  Results <- append(Results, list("ValErrs0" = ValErrs0))
}
Results <- append(Results, compResults(MCIS0, SigmaArray, "MCIS0"))



# ------------------------------------
# Oracle
# ------------------------------------
oracleSigmaSoft <- array(0, dim=c(J,p,p))
oracleSigmaHard <- array(0, dim=c(J,p,p))
for(j in 1:J){
	oracleSigmaSoft[j,,] <- adaptThresoldCov(Y[j,,], Yval[j,,], soft = 1)$sigma
	oracleSigmaHard[j,,] <- adaptThresoldCov(Y[j,,], Yval[j,,], soft = 0)$sigma
	
}
Results <- append(Results, compResults(oracleSigmaSoft, SigmaArray, "OracleSoft"))
Results <- append(Results, compResults(oracleSigmaHard, SigmaArray, "OracleHard"))



#-------------------------------------
# cclasso-val
# ------------------------------------
cclassoEst <- array(0, dim=c(J,p,p))
for(j in 1:J){
	t1 <- cclasso(X[j,,], counts = FALSE, pseudo = 0.5, k_cv = 10, 
                    lam_int = c(1e-4, 1), k_max = 500, n_boot = 50)
	cclassoEst[j,,] <- diag(sqrt(t1$var_w))%*%t1$cor_w%*%diag(sqrt(t1$var_w))
}

Results <- append(Results, compResults(cclassoEst, SigmaArray, "cclasso"))
saveRDS(Results, savename)

rm(list = ls())
q("no")
