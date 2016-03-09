#######################################################################
# Airport Clustering Research                                         #
# Data was attained from RITA for all available airports in July 2015 #
# I can provide the database on request, evangfreitag-at-gmail.com    #
#######################################################################
# For plotting
require(ggplot2)
require(seriation)
require(plotly)
require(MASS)
require(scatterplot3d)

require(cluster) # For Gap statistics
require(LICORS) # For kmeans++

############################
# Step 1: Data Preparation #
############################
	###############################
	# A. Open the SQLite database #
	###############################
	setwd("C:/Users/Evan/Documents/AirlineResearch/Clustering")
	require(RSQLite)
	Data <- dbConnect(SQLite(), dbname = "AirlineResearch.db")
	initExtension(Data) # Loads special math functions such the as median
	from_db <- function(sql) {dbGetQuery(Data, sql)}

	#######################
	# B. Read in the data #
	#######################
	# The data that is going to be analyzed
	Delay <- from_db("select DAY_OF_MONTH, UNIQUE_CARRIER, ORIGIN, 
				DEP_TIME, WEATHER_DELAY, NAS_DELAY, SECURITY_DELAY, CARRIER_DELAY,
				LATE_AIRCRAFT_DELAY, DEP_DELAY_MINUTES, airport
				from AR2015 
				inner join airports on AR2015.ORIGIN = airports.iata and DEP_DELAY_MINUTES > 0
			     	order by ORIGIN, DAY_OF_MONTH, DEP_TIME")


	AirportData <- from_db("select distinct ORIGIN, airport
				from AR2015 
				inner join airports on AR2015.ORIGIN = airports.iata and DEP_DELAY_MINUTES > 0
			     	order by ORIGIN")
		
	rm(Data, from_db)

	#####################################################################
	# The missing values in the data were replaced with -9999 initially #
	# Replacing the -9999's with 0's as a blank for any type of delay   #
	# time seems to indicate that the delay was not that type           #                                    #
	#####################################################################
	Delay[Delay == -9999] <- 0

	#######################################
	# C. Split the data by origin airport #
	#######################################
	# Note: Output is a list where each element is the data for a unique origin
	SPLIT_BY_ORIGIN <- split(Delay, Delay$ORIGIN)
		
	##################################################################
	# D. Pick out the vectors of delay times and save them in a list #
	##################################################################
	DELAY_BY_ORIGIN_LIST <- vector("list", length(SPLIT_BY_ORIGIN))
	for(n in 1:length(SPLIT_BY_ORIGIN)){
		DELAY_BY_ORIGIN_LIST[[n]] <- SPLIT_BY_ORIGIN[[n]][, c("DAY_OF_MONTH", 
										"WEATHER_DELAY", 
										"NAS_DELAY", 
										"SECURITY_DELAY", 
										"CARRIER_DELAY", 
										"LATE_AIRCRAFT_DELAY")]
	}

############################################
# Step 2: Calculate Dissimilarity Matrices #
############################################
	###############################
	# Multivariate CID SVD MaxMin #
	###############################
	# MCIDSVDMaxMin <- MultiDissimilaritySVDMaxMin(DELAY_BY_ORIGIN_LIST, DistMethod = "CID")

	###############################
	# Multivariate DTW SVD MaxMin #
	###############################
	# MDTWSVDMaxMin <- MultiDissimilaritySVDMaxMin(DELAY_BY_ORIGIN_LIST, DistMethod = "DTW")

	##################################
	# Multivariate CIDDTW SVD MaxMin #
	##################################
	# MCIDDTWSVDMaxMin <- MultiDissimilaritySVDMaxMin(DELAY_BY_ORIGIN_LIST, DistMethod = "CIDDTW")

	############################################
	# Name the DELAY_BY_ORIGIN_LIST components #
	############################################
	# This happens after the dissimilarity computations 
	# as named lists are treated differently in R.
	# Non-named lists are faster to process.
	names(DELAY_BY_ORIGIN_LIST) <- c(names(SPLIT_BY_ORIGIN))
	rm(SPLIT_BY_ORIGIN)

	#################################################################################
	# Define the addNames function to pass list names to the dissimilarity matrices #
	#################################################################################
	# addNames <- function(matrix, list) {
	#	colnames(matrix) <- names(list)
	#	rownames(matrix) <- names(list)
	#	return(matrix)
	#}

	#############################################################
	# Add the names to the dissimilarity matrices and save them #
	#############################################################
	# MCIDSVDMaxMin <- addNames(MCIDSVDMaxMin, DELAY_BY_ORIGIN_LIST)
	# MDTWSVDMaxMin <- addNames(MDTWSVDMaxMin, DELAY_BY_ORIGIN_LIST)
	# MCIDDTWSVDMaxMin <- addNames(MCIDDTWSVDMaxMin, DELAY_BY_ORIGIN_LIST)
	# save(MCIDSVDMaxMin, file = "MCIDSVDMaxMin.RData")
	# save(MDTWSVDMaxMin, file = "MDTWSVDMaxMin.RData")
	# save(MCIDDTWSVDMaxMin, file = "MCIDDTWSVDMaxMin.RData")

	#########################################
	# Load the saved Dissimilarity Matrices #
	#########################################
	load("MCIDSVDMaxMin.RData")
	load("MDTWSVDMaxMin.RData")
	load("MCIDDTWSVDMaxMin.RData")

###########################################
# Step 3: Plot the Dissimilarity Matrices #
###########################################
	#######################################################################
	# Various plots for Multivariate CID SVD Max-Min dissimilarity matrix #
	#######################################################################
	PlotDissimAsVector(MCIDSVDMaxMin, Title=c("Ordered MCIDSVD Max-Min Dissimilarity Matrix","MCIDSVD Max-Min Dissimilarity Matrix"))
	RelativePlots(MCIDSVDMaxMin, Names=names(ORIGIN_LIST), Ylabel="MCIDSVD Max-Min Dissimilarity from All Other Origins", Title ="Relative MCIDSVD Max-Min Dissimilarities")
	ColorDissimilarityPlots(MCIDSVDMaxMin, maxColors=10, byrank=TRUE, Title= c("MCIDSVD Max-Min Dissimilarity Matrix","MCIDSVD Max-Min Dissimilarity Matrix"))
	ColorDissimilarityPlots(MCIDSVDMaxMin, maxColors=10, byrank=FALSE, Title=c("MCIDSVD Max-Min Dissimilarity Matrix","MCIDSVD Max-Min Dissimilarity Matrix"))
	SeriationPlots(MCIDSVDMaxMin, "MCIDSVD Max-Min Dissimilarity Matrix", maxK=10)
	HeatmapPlot(MCIDSVDMaxMin)
	ScatterMDS2D(MCIDSVDMaxMin)
	ScatterMDS3D(MCIDSVDMaxMin)
	SurfaceMDS3D(MCIDSVDMaxMin)

	#######################################################################
	# Various plots for Multivariate DTW SVD Max-Min dissimilarity matrix #
	#######################################################################
	PlotDissimAsVector(MDTWSVDMaxMin, Title=c("Ordered MDTWSVD Max-Min Dissimilarity Matrix","MDTWSVD Max-Min Dissimilarity Matrix"))
	RelativePlots(MDTWSVDMaxMin, Names=names(ORIGIN_LIST), Ylabel="MDTWSVD Max-Min Dissimilarity from All Other Origins", Title ="Relative MDTWSVD Max-Min Dissimilarities")
	ColorDissimilarityPlots(MDTWSVDMaxMin, maxColors=10, byrank=TRUE, Title= c("MDTWSVD Max-Min Dissimilarity Matrix","MDTWSVD Max-Min Dissimilarity Matrix"))
	ColorDissimilarityPlots(MDTWSVDMaxMin, maxColors=10, byrank=FALSE, Title=c("MDTWSVD Max-Min Dissimilarity Matrix","MDTWSVD Max-Min Dissimilarity Matrix"))
	SeriationPlots(MDTWSVDMaxMin, "MDTWSVD Max-Min Dissimilarity Matrix", maxK=10)
	HeatmapPlot(MDTWSVDMaxMin)
	ScatterMDS2D(MDTWSVDMaxMin)
	ScatterMDS3D(MDTWSVDMaxMin)
	SurfaceMDS3D(MDTWSVDMaxMin)

	##########################################################################
	# Various plots for Multivariate CIDDTW SVD Max-Min dissimilarity matrix #
	##########################################################################
	PlotDissimAsVector(MCIDDTWSVDMaxMin, Title=c("Ordered MCIDDTWSVD Max-Min Dissimilarity Matrix","MCIDDTWSVD Max-Min Dissimilarity Matrix"))
	RelativePlots(MCIDDTWSVDMaxMin, Names=names(ORIGIN_LIST), Ylabel="MCIDDTWSVD Max-Min Dissimilarity from All Other Origins", Title ="Relative MCIDDTWSVD Max-Min Dissimilarities")
	ColorDissimilarityPlots(MCIDDTWSVDMaxMin, maxColors=10, byrank=TRUE, Title= c("MCIDDTWSVD Max-Min Dissimilarity Matrix","MCIDDTWSVD Max-Min Dissimilarity Matrix"))
	ColorDissimilarityPlots(MCIDDTWSVDMaxMin, maxColors=10, byrank=FALSE, Title=c("MCIDDTWSVD Max-Min Dissimilarity Matrix","MCIDDTWSVD Max-Min Dissimilarity Matrix"))
	SeriationPlots(MCIDDTWSVDMaxMin, "MCIDDTWSVD Max-Min Dissimilarity Matrix", maxK=10)
	HeatmapPlot(MCIDDTWSVDMaxMin)
	ScatterMDS2D(MCIDDTWSVDMaxMin)
	ScatterMDS3D(MCIDDTWSVDMaxMin)
	SurfaceMDS3D(MCIDDTWSVDMaxMin)

####################################
# Step 4: Calculate Gap Statistics #
####################################
# 1,000 bootstrap iterations
# Range of k is from 1 to 50 (1 as an optimum indicates no clusters)

	##################################
	# Function to get various optima #
	##################################
	GapValues <- function(x){
		k1 <- maxSE(x$Tab[, "gap"], x$Tab[, "SE.sim"], method="Tibs2001SEmax")
		k2 <- maxSE(x$Tab[, "gap"], x$Tab[, "SE.sim"], method="firstSEmax")
		k3 <- maxSE(x$Tab[, "gap"], x$Tab[, "SE.sim"], method="globalSEmax")
		k4 <- maxSE(x$Tab[, "gap"], x$Tab[, "SE.sim"], method="firstmax")
		k5 <- maxSE(x$Tab[, "gap"], x$Tab[, "SE.sim"], method="globalmax")
	return(list(k1,k2,k3,k4,k5))
	}

	####################################
	# Gap Statistics for MCIDSVDMaxMin #
	####################################
	gap_MCIDSVDMaxMin <- clusGap(MCIDSVDMaxMin, FUN = kmeans, nstart = 1, K.max=50, B = 1000, verbose = interactive())
	plot(gap_MCIDSVDMaxMin)
	k_MCIDSVDMaxMin <- GapValues(gap_MCIDSVDMaxMin)
	k_MCIDSVDMaxMin
	
	####################################
	# Gap Statistics for MDTWSVDMaxMin #
	####################################
	gap_MDTWSVDMaxMin <- clusGap(MDTWSVDMaxMin, FUN = kmeans, nstart = 1, K.max=50, B = 1000, verbose = interactive())
	plot(gap_MDTWSVDMaxMin)
	k_MDTWSVDMaxMin <- GapValues(gap_MDTWSVDMaxMin)
	k_MDTWSVDMaxMin
	
	#######################################
	# Gap Statistics for MCIDDTWSVDMaxMin #
	#######################################
	gap_MCIDDTWSVDMaxMin <- clusGap(MCIDDTWSVDMaxMin, FUN = kmeans, nstart = 1, K.max=50, B = 1000, verbose = interactive())
	plot(gap_MCIDDTWSVDMaxMin)
	k_MCIDDTWSVDMaxMin <- GapValues(gap_MCIDDTWSVDMaxMin)
	k_MCIDDTWSVDMaxMin

#######################
# Calculate K-Means++ #
#######################
	###############################
	# K-Means++ for MCIDSVDMaxMin #
	###############################
	km_MCIDSVDMaxMin <- kmeanspp(MCIDSVDMaxMin, 5, iter.max=500)
	ClusterNumberCID <- km_MCIDSVDMaxMin$cluster
	
	###############################
	# K-Means++ for MDTWSVDMaxMin #
	###############################
	km_MDTWSVDMaxMin <- kmeanspp(MDTWSVDMaxMin, 5, iter.max=500)
	ClusterNumberDTW <- km_MDTWSVDMaxMin$cluster

	##################################
	# K-Means++ for MCIDDTWSVDMaxMin #
	##################################
	km_MCIDDTWSVDMaxMin <- kmeanspp(MCIDDTWSVDMaxMin, 16, iter.max=500)
	ClusterNumberCIDDTW <- km_MCIDDTWSVDMaxMin$cluster

	#########################################
	# Bind the clusters to the airport data #
	#########################################
	FinalClusters <- cbind(AirportData, ClusterNumberCID, ClusterNumberDTW, ClusterNumberCIDDTW)

	##############################################
	# Merge the Delay, airport, and cluster data #
	##############################################
	FinalClusters <- merge(Delay, FinalClusters, by.Delay = "ORIGIN", by.FinalClusters = "ORIGIN")
	FinalClustersList <- split(FinalClusters, FinalClusters$ClusterNumberCIDDTW)

	###################################################################
	# Function to get the airport names contained within each cluster #
	###################################################################
	GetNames <- function(alist) {
		InClusterNames <- list()
		for (i in 1:length(alist)) {
			InClusterNames[[i]] <- unique(alist[[i]]$airport)
		}
	return(InClusterNames)
	}

	#######################################################
	# Get the airport names contained within each cluster #
	#######################################################
	Names <- GetNames(FinalClustersList)

	###############################################################
	# Function to split the cluster lists by origin and name them #
	###############################################################
	GetLofL <- function(alist, Names) {
		ListOfClusterLists <- list()
		for (i in 1:length(alist)) {
			ListOfClusterLists[[i]] <- split(alist[[i]], alist[[i]]$ORIGIN)
			names(ListOfClusterLists[[i]]) <- Names[[i]]
		}
	return(ListOfClusterLists)
	}

	###################################################
	# Split the cluster lists by origin and name them #
	###################################################
	FinalClustersList <- GetLofL(FinalClustersList, Names)

	###############################################################################################################
	# Remove the clusterings and columns we're not interested in, in this case, only looking at CIDDTW clustering #
	###############################################################################################################
	for(i in 1:length(FinalClustersList)){
		for(j in 1:length(FinalClustersList[[i]])){
			FinalClustersList[[i]][[j]]$airport <- NULL
			FinalClustersList[[i]][[j]]$ORIGIN <- NULL
			FinalClustersList[[i]][[j]]$UNIQUE_CARRIER <- NULL
			FinalClustersList[[i]][[j]]$ClusterNumberCID <- NULL
			FinalClustersList[[i]][[j]]$ClusterNumberDTW <- NULL
			FinalClustersList[[i]][[j]]$DEP_TIME <- NULL
		}
	}

#####################
# Plot the clusters #
#####################
# To visualize the clusters, select a cluster and plot with the MultivariateTimeSeriesPlots function

# Define colors as this prints all time series for each individual airport
Colors <- c("blue", "green", "orange", "yellow", "red", "black", "gray")

# LegendNames <- c("Weather Delay", "NAS Delay", "Security Delay", "Carrier Delay", "Late Aircraft Delay", "Departure Delay", "No Responsibility")
MultivariateTimeSeriesPlots(FinalClustersList[[1]], MarkerType="b", Xlabel="Day of Month", Xlimit=c(1, 31), Ylimit="TotalMaximum", PointSize=.75, LineWidth=2, Names=names(ORIGIN_LIST), Colors)
MultivariateTimeSeriesPlots(FinalClustersList[[2]], MarkerType="b", Xlabel="Day of Month", Xlimit=c(1, 31), Ylimit="TotalMaximum", PointSize=.75, LineWidth=2, Names=names(ORIGIN_LIST), Colors)
MultivariateTimeSeriesPlots(FinalClustersList[[3]], MarkerType="b", Xlabel="Day of Month", Xlimit=c(1, 31), Ylimit="TotalMaximum", PointSize=.75, LineWidth=2, Names=names(ORIGIN_LIST), Colors)
MultivariateTimeSeriesPlots(FinalClustersList[[4]], MarkerType="b", Xlabel="Day of Month", Xlimit=c(1, 31), Ylimit="TotalMaximum", PointSize=.75, LineWidth=2, Names=names(ORIGIN_LIST), Colors)
MultivariateTimeSeriesPlots(FinalClustersList[[5]], MarkerType="b", Xlabel="Day of Month", Xlimit=c(1, 31), Ylimit="TotalMaximum", PointSize=.75, LineWidth=2, Names=names(ORIGIN_LIST), Colors)

