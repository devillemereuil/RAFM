l5 <-
function(kap_, z_, dat,dom_,Fis_){
	K = matrix(kap_, nrow=nrow(z_), ncol=ncol(z_))
	p = apply(K*z_, 2, sum)
	if (!dom_)	{
	  f = dmultinom(dat, size=sum(dat), prob=p, log=TRUE)
	} else { #From Foll & Gaggiotti (2008)
	  g2=(1-Fis_)*(1-p[2])**2 + Fis_*(1-p[2])
	  f = dbinom(dat[2], size=sum(dat), prob=1-g2, log=TRUE)
	}
 	return(f)
	}

