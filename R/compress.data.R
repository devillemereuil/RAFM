compress.data <-
function(rawDat,dom=FALSE,unobs=FALSE){
	pops = as.factor(rawDat[,1])
	DNA = rawDat[,2:ncol(rawDat)]
	if (dom) unobs=FALSE
	if (!dom) {
	  dna = matrix(NA, ncol=ncol(DNA)/2, nrow=2*nrow(DNA))
	  for( i in 1:(ncol(DNA)/2) ){
		  dna[,i] = c(DNA[,2*i], DNA[,-1+2*i])
		}
	} else {
	  dna<-DNA
	}
	pops = rep(pops, 2)
	uniqs = levels(pops)
	nloci = ncol(dna)
	dat = list(1:nloci)
	for( j in 1:nloci ){
		dnaj = dna[,j]
		nas = which(is.na(dnaj))
		if( length(nas) > 0 ){
			dnaj = dnaj[-nas]
			popsj = pops[-nas]
			}
		else{ popsj = pops }
		if (unobs) {genotypes = c(unique(dnaj), "unobs")} else {genotypes = unique(dnaj)}
		datj = matrix(0, nrow=length(uniqs), ncol=length(genotypes))
    rownames(datj) = uniqs
		for( p in 1:length(uniqs) ){
			thispop = dnaj[which(popsj==uniqs[p])]
			datj[uniqs[p], ] = sapply(genotypes, function(x) length(which(thispop==x)))
			}
		dat[[j]] = datj
		}
	return(dat)
	}

