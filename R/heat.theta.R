heat.theta <-
function(thetapost, mean=TRUE, verbose=TRUE){
	mod = as.matrix(0*thetapost[,,1])
	central = lower = upper = mod
	npop = dim(thetapost)[1]
	for( i in 1:npop ){
		for( j in 1:npop ){
			X = thetapost[i,j,]
			lower[i,j] = lower[j,i] = quantile(X, 0.025)
			upper[i,j] = upper[j,i] = quantile(X, 0.975)
			if( mean ){ 
				central[i,j] = central[j,i] = mean(X)
				}
			else{ central[i,j] = central[j,i] = median(X) }
			}
		}
	central = 0.5*(central + t(central))
	lower = 0.5*(lower + t(lower))
	upper = 0.5*(upper + t(upper))
	out = as.list(1:3)
	out[[1]] = central
	out[[2]] = lower
	out[[3]] = upper
	names(out) = c("mean", "lower", "upper")
	if( !mean ){
		names(out)[1] = "median"
		}
	heatmap(central, symm=T, margins=c(5,5), xlab="subpopulations, note permutation", ylab="", main="Coancestry matrix")
	if( verbose ){
		return(out)
		}
	}

