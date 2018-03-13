do.all <-
function(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1, dom=FALSE,unobs=FALSE,plot=FALSE){
  primary <- AFM(dat=dat, nMC=nMC, burnin=burnin, thinning=thinning, eps=eps, priq=priq, pria=pria, prik=prik, dom=dom,unobs=unobs)
	kap <- primary[[1]]
	alpha <- primary[[2]]
	secondary <- gen.theta(kap, alpha)
	theta <- secondary[[1]]
	fst <- secondary[[2]]
	alpha <- exp(alpha)
	npop = dim(theta)[1]
  if (plot) {
  	print("Drawing trace of coancestry matrix theta, assess convergence.", quote=F)
  	if( npop > 4 ){ 
  		npop = 4
  		print("N.B. Omitting populations 5+.", quote=F)
  		}
  	# windows()
  	par(mfcol=c(npop,npop))
  	for( i in 1:npop ){
  		for( j in 1:npop ){
  			headl = paste("theta[ ",i," , ",j," ]", sep="")
  			plot( theta[i,j,], type='l', main=headl, xlab="", ylab="" )
  			}
    }
  }
	out = as.list(1:4)
	names(out) = c("theta", "fst", "kappa", "alpha")
	out[[1]] = theta
	out[[2]] = fst
	out[[3]] = kap
	out[[4]] = alpha
	return(out)
	}
