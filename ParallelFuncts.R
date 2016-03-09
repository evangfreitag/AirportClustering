#############################################################
# Functions to calculate Dissimilarity Measures in parallel #
#############################################################
require(foreach)
n_cores <- detectCores() -1
cl <- makeCluster(n_cores)

###########################################################
# Define the Complexity Invariant Distance (CID) function #
###########################################################
cid <- function(x,y) {
	cesx <- sqrt(sum(diff(x)^2))
	cesy <- sqrt(sum(diff(y)^2))
	denom <- min(cesx, cesy)
	if(denom == 0){denom <- 1}
	cid1 <- max(cesx, cesy)/denom 
	return(cid1)
}

########################################	
# Define the Multivariate CID function #
########################################
CidMat <- function(x,y) {
	N <- ncol(x)
	m1 <- matrix(0,N,N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:N) %dopar% {
			m1[i,j] <- cid(x[,i], y[,j])
	}
	return(m1)
}
	
########################################
# Define the Multivariate DTW function #
########################################
DtwMat <- function(x,y) {
	require(dtw)
	N <- ncol(x)
	m1 <- matrix(0,N,N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:N) %dopar% {
			m1[i,j] <- dtw(x[,i], y[,j], distance.only = "TRUE")$distance
	}
	return(m1)	
}

###################################################
# Define the Normalized Multivariate DTW function #
###################################################
NormalizedDtwMat <- function(x,y) {
	require(dtw)
	N <- ncol(x)
	m1 <- matrix(0,N,N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:N) %dopar% {
			m1[i,j] <- dtw(x[,i], y[,j], distance.only = "TRUE")$normalized
	}
	return(m1)
}	

###################################################
# Define the Multivariate DTW SVD MaxMin function #
###################################################
Multi_DTW_SVD_MaxMin <- function(alist) {
	# Get the list length
	N <- length(alist)
	m1 <- matrix(0, N, N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:(i-1)) %dopar% {
			m1[i,j] <- svd(DtwMat(alist[[i]],alist[[j]]))$d
	}
	return(m1 + t(m1))
}

###################################################
# Define the Multivariate CID SVD MaxMin function #
###################################################
Multi_CID_SVD_MaxMin <- function(alist) {
	# Get the list length
	N <- length(alist)
	m1 <- matrix(0, N, N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:(i-1)) %dopar% {
			m1[i,j] <- svd(CidMat(alist[[i]],alist[[j]]))$d
	}
	return(m1 + t(m1))
}
	
##############################################################
# Define the Multivariate Normalized DTW SVD MaxMin function #
##############################################################	
Multi_NormalizedDTW_SVD_Max_Min <- function(alist) {
	# Get the list length
	N <- length(alist)
	m1 <- matrix(0, N, N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:(i-1)) %dopar% {
			m1[i,j] <- svd(NormalizedDtwMat(alist[[i]],alist[[j]]))$d
	}
	return(m1 + t(m1))
}

######################################################
# Define the Multivariate CIDDTW SVD MaxMin function #
######################################################
Multi_CIDDTW_SVD_MaxMin <- function(alist) {
	# Get the list length
	N <- length(alist)
	m1 <- matrix(0, N, N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:(i-1)) %dopar% {
			m1[i,j] <- svd(DtwMat(alist[[i]],alist[[j]])*CidMat(alist[[i]],alist[[j]]))$d
		}
	return(m1 + t(m1))
}

##################################################################
# Define the Multivariate CID Normalized DTW SVD MaxMin function #
##################################################################
Multi_CIDNormalizedDTW_SVD_MaxMin <- function(alist) {
	# Get the list length
	N <- length(alist)
	m1 <- matrix(0, N, N)
	x1 <- foreach(i=1:N) %:%
		foreach(j=1:(i-1)) %dopar% {
			m1[i,j] <- svd(NormalizedDtwMat(alist[[i]],alist[[j]])*CidMat(alist[[i]],alist[[j]]))$d
		}
	return(m1 + t(m1))
}

