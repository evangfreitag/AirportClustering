###################################################################################
# Function to calculate a dissimilarity matrix for multivariate data based on svd #
###################################################################################
# Step 1 is to calculate individual dissimilarity matrices for each pair of matrices
# Step 2 is to calculate the svd of the individual dissimilarity matrices
# Step 3 is to map the d in the svd from R to R^n using max-min for each dissimilarity matrix
# Step 4 is to store the mapped values in a new (single) dissimilarity matrix 
# The input is a list of matrices

MultiDissimilaritySVDMaxMin <- function(alist, DistMethod) {

	if(DistMethod == "CID" | DistMethod == "CIDDTW"){
		# Define the Complexity Invariant Distance (CID) function
		cid <- function(x,y) {
			cesx <- sqrt(sum(diff(x)^2))
			cesy <- sqrt(sum(diff(y)^2))
			denom <- min(cesx, cesy)
			if(denom == 0){denom <- 1}
			cid1 <- max(cesx, cesy)/denom 
			return(cid1)
			}
	
		# Define the CidMat function
		CidMat <- function(matrix1,matrix2) {
			X <- as.matrix(matrix1)
			Y <- as.matrix(matrix2)
			N <- length(X[1,])
			CidM <- rep(0, N)
			CidM2 <- matrix(0, N, N)
			for(i in 1:N){
				for(j in 1:N){
					CidM[j] <- cid(X[,i],Y[,j])
				}
				CidM2[i,] <- CidM
			}
			return(as.matrix(CidM2))
		}
	}

	if(DistMethod == "DTW" | DistMethod == "CIDDTW"){
		require(dtw)
		# Define the DtwMat function
		DtwMat <- function(matrix1, matrix2) {
			X <- as.matrix(matrix1)
			Y <- as.matrix(matrix2)
			N <- length(X[1,])
			DtwM <- rep(0, N)
			DtwM2 <- matrix(0, N, N)
			for(i in 1:N){
				for(j in 1:N){
					DtwM[j] <- dtw(X[,i],Y[,j], distance.only = TRUE)$distance
				}
			DtwM2[i,] <- DtwM
			}
		return(as.matrix(DtwM2))
		}
	}

	if(DistMethod == "NormalizedDTW" | DistMethod == "CIDNormalizedDTW"){
		require(dtw)
		# Define the DtwMat function
		NormalizedDtwMat <- function(matrix1, matrix2) {
			X <- as.matrix(matrix1)
			Y <- as.matrix(matrix2)
			N <- length(X[1,])
			DtwM <- rep(0, N)
			DtwM2 <- matrix(0, N, N)
			for(i in 1:N){
				for(j in 1:N){
					DtwM[j] <- dtw(X[,i],Y[,j], distance.only = TRUE)$normalized
				}
			DtwM2[i,] <- DtwM
			}
		return(as.matrix(DtwM2))
		}
	}
	
	# Get the list length
	N <- length(alist)
	
	# Preallocate the outer loop list to hold the results from the inner loop
	OuterMatrix <- matrix(0, N, N)
	InnerVector <- rep(0, N)
	
		# Start the outer part of the nested loop
		for (i in 1:N) {

			# Preallocate the list and vector within the inner loop
			InnerList <- vector("list", i)	

			# Preallocate the matrices within the inner lists
			for(k in 1:i){
				InnerList[[k]] <- matrix(0,N,N)
			}			

				# Start the inner loop				
				for (j in 1:i) {
					if(j < i){
						# Get the dissimilarity matrix for the jth list member in regards to the ith
						if(DistMethod == "CID") {InnerList[[j]] <- svd(CidMat(alist[[i]],alist[[j]]))$d}
						else if(DistMethod == "DTW") {InnerList[[j]] <- svd(DtwMat(alist[[i]],alist[[j]]))$d}
						else if(DistMethod == "NormalizedDTW") {InnerList[[j]] <- svd(NormalizedDtwMat(alist[[i]],alist[[j]]))$d}
						else if(DistMethod == "CIDDTW") {InnerList[[j]] <- svd(CidMat(alist[[i]],alist[[j]])*DtwMat(alist[[i]],alist[[j]]))$d}
						else if(DistMethod == "CIDNormalizedDTW") {InnerList[[j]] <- svd(CidMat(alist[[i]],alist[[j]])*NormalizedDtwMat(alist[[i]],alist[[j]]))$d}

						# Reduce the d in svd
						InnerVector[j] <- max(InnerList[[j]])-min(InnerList[[j]])
					}					
				}

		# Save the output from the inner part of the nested loop
		OuterMatrix[,i] <- InnerVector
		
		# Print the progress in terms of i
		print(i)
		}
	
	# Unlist the outer list
	DissimilarityMatrix <- OuterMatrix + t(OuterMatrix)

	return(DissimilarityMatrix)
}