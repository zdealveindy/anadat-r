'PCNM' <- 
function(matdist, thresh=NULL, dbMEM=FALSE, moran=NULL, all=FALSE, include.zero=FALSE, silent=FALSE)
#
# Compute the PCNM or dbMEM eigenfunctions corresponding to 
# all eigenvalues (+, 0, -). 
#    In PCNM computation, the diagonal of D = 0.
#    In dbMEM, the diagonal of D = 4*threshh.
#    Distance-based MEM are described in Dray et al. 2006. 
#    The name was abbreviated to db-MEM by PPN & PL (subm.)
# Input file: distance matrix produced by the function "dist".
# Computation of the threshold requires a function of the library "ape".
# 
# Original PCNM function: Stephane Dray, November 11, 2004
# The present version: Pierre Legendre, August 2007, January and March 2009
{
	require(vegan)
	epsilon <- sqrt(.Machine$double.eps) 
	a <- system.time({
	if(is.null(moran)) {
		if(dbMEM) { moran=FALSE } else { moran=TRUE }
		}
	single <- FALSE
	if(moran) {
		# cat("The site coordinates were computed from 'matdist'.",'\n')
		pcoa.xy <- pcoa.all(matdist)
		
		if(is.na(pcoa.xy$values[2]) | (pcoa.xy$values[2] < epsilon)) {
			if(!silent) cat("The sites form a straight line on the map.",'\n')
			xy <- pcoa.xy$vectors
			single <- TRUE
			} else {
			xy <- pcoa.xy$vectors[,1:2]
			}
		}

	matdist <- as.matrix(matdist)
	n <- nrow(matdist)

	# Truncation of distance matrix
	if(is.null(thresh)) {
		spanning <- vegan::spantree(as.dist(matdist))
		threshh <- max(spanning$dist)
	    if(!silent) cat("Truncation level =",threshh+0.000001,'\n')
		} else {
		threshh = thresh
	    if(!silent) cat("User-provided truncation threshold =",thresh,'\n')
		}
	matdist[matdist > threshh] <- 4*threshh

	if(dbMEM==FALSE) { diagonal <- 0 } else { diagonal <- 4*threshh }

	mypcnm.all <- pcoa.all(matdist, diagonal=diagonal, all=all, include.zero=include.zero, rn=rownames(matdist))

	# Compute Moran's I
	if(moran) {
		require(AEM)
		if(single) {
			nb <- dnearneigh(matrix(c(xy,rep(0,n)),n,2), 0, (threshh + epsilon))
			} else {
			nb <- dnearneigh(xy, 0, (threshh + epsilon))
			}
		fr.to.pcnm2 <- as.matrix(listw2sn(nb2listw(nb))[,1:2])
		weight.dist.coord.mat <- as.matrix(1-(as.dist(matdist)/(4*threshh))^2)
		weight <- weight.dist.coord.mat[fr.to.pcnm2]
		res <- moran.I.multi(mypcnm.all$vectors, link=fr.to.pcnm2, weight=weight)
		Moran <- res$res.mat[,1:2]
		positive <- rep(FALSE,length(mypcnm.all$values))
		positive[which(Moran[,1] > res$expected)] <- TRUE
		Moran <- cbind(as.data.frame(Moran), positive)
		colnames(Moran) <- c("Moran","p.value","Positive")
		}
	})
	a[3] <- sprintf("%2f",a[3])
	if(!silent) cat("Time to compute PCNMs =",a[3]," sec",'\n')
	if(is.null(thresh)) {
		if(moran) {
			res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, spanning=spanning, thresh=threshh+0.000001)
			} else {
			res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, spanning=spanning, thresh=threshh+0.000001)
			}
		} else {
		if(moran) {
			res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, thresh=thresh)
			} else {
			res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, thresh=threshh+0.000001)
			}
		}
	res
}
