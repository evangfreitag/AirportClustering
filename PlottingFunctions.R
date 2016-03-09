######################
# Plotting Functions #
######################
###############################################################
# The readkey function: pauses a for loop until you hit enter #
###############################################################
readkey <- function() {
   	cat ("Press [enter] to continue")
	line <- readline()
}

################################################################################
# MultivariateTimeSeriesPlots: Plot each list element individually with the    #
# y-axis constant (YLIM <-"Total") or with the y-axis dependent on the current #
# list element (YLIM="Individual")	                                       #
################################################################################
MultivariateTimeSeriesPlots <- function(alist, MarkerType, Xlabel, Xlimit, Ylimit, PointSize, LineWidth, Names, Colors){
	if(Ylimit == "TotalMaximum"){ YMax <- max(sapply(alist, max))
					   YLIM <- c(0, YMax)}
	for(i in 1:length(alist)) {
		if(Ylimit == "IndividualMaximum"){YLIM <- c(0,max(alist[[i]][,-1]))}
		Name <- names(alist)[i]
		for(j in 2:length(Colors)){
			plot(alist[[i]][,1], alist[[i]][,j], type=MarkerType, lwd=LineWidth, cex=PointSize, col=Colors[j], xlim=Xlimit, ylim=YLIM, xlab = Xlabel, ylab = "", main = Name)
			par(new=TRUE)
		}
	readkey()
	par(new=FALSE)
	}
}

######################################################################################
# PlotDissimAsVectors: Plot a list of dissimilarity matrices individually as vectors #
# Vectors can be sorted and/or scaled if desired                                     #
######################################################################################
PlotDissimAsVectors <- function(DissimilaritiesList, SORT, SCALE, Titles){
	if(SCALE == "FALSE"){ 
		YMax <- max(sapply(DissimilaritiesList, max))
		YLIM <- c(0, YMax)
	}	
	if(SORT == "TRUE" & SCALE == "TRUE"){
		for(i in 1:length(Titles)) {
			plot(sort(scale(as.numeric(DissimilaritiesList[[i]]))), cex=.25, ylab="Dissimilarity Values", main=Titles[i])
			par(new=TRUE)
		}
	}
	if(SORT == "TRUE" & SCALE == "FALSE"){
		for(i in 1:length(Titles)) {
			plot(sort(as.numeric(DissimilaritiesList[[i]])), cex=.25, ylim = YLIM, ylab="Dissimilarity Values", main=Titles[i])
			par(new=TRUE)
		}
	}
	if(SORT == "FALSE" & SCALE == "TRUE"){
		for(i in 1:length(Titles)) {
			plot(scale(as.numeric(DissimilaritiesList[[i]])), cex=.25, ylab="Dissimilarity Values", main=Titles[i])
			par(new=TRUE)
		}
	}
	if(SORT == "FALSE" & SCALE == "FALSE"){
		for(i in 1:length(Titles)) {
			plot(as.numeric(DissimilaritiesList[[i]]), cex=.25, ylim = YLIM, ylab="Dissimilarity Values", main=Titles[i])
			par(new=TRUE)
		}
	}
}

###############################################################
# PlotDissimAsVector: Plot a dissimilarity matrix as a vector #
# plots the vector both sorted and unsorted                   # 
###############################################################
PlotDissimAsVector <- function(DissimilarityMatrix, Titles){
	par(mfrow=c(1,2))
	plot(sort(as.numeric(DissimilarityMatrix)), cex=.25, ylab="Dissimilarity Values", xlab="", main=Titles[1])
	plot(as.numeric(DissimilarityMatrix), cex=.25, ylab="Dissimilarity Values", xlab="", main=Titles[2])
}

#############################################################################################
# RelativePlots: Plot each element in relation to all other elements based on dissimilarity #
#############################################################################################
RelativePlots <- function(DissimilarityMatrix, Names, Ylabel, Title){
	Ylim <- c(0, ceiling(max(DissimilarityMatrix)))
	for(i in 1:length(DissimilarityMatrix[,1])){
		Colors <- rep("orange", length(DissimilarityMatrix[,1]))
		Colors[i] <- "blue"
		plot(DissimilarityMatrix[,i], ylim=Ylim, xlab=Names[i], ylab=Ylabel, col=Colors, pch=19, main=Title)
	readkey()
	}
}

###################################################################################
# ColorDissimilarityPlot: Plots for unordered and ordered dissimilarity matrices  #
# maxColors sets number of colors/clusters                                        #
###################################################################################
ColorDissimilarityPlots <- function(DissimilarityMatrix, maxColors, byrank, Title){
	
	# Color dissimilarity function
	ColorDissimilarity <- function(DissimilarityMatrix, NumberOfColors, byrank, Title){
		require(gclus)
		if(max(DissimilarityMatrix) > 1) DissimilarityMatrix <- DissimilarityMatrix/max(DissimilarityMatrix)
		if(byrank == TRUE) {colors1 = dmat.color(1-DissimilarityMatrix, cm.colors(NumberOfColors))}		
		else {colors1 = dmat.color(1-DissimilarityMatrix, byrank=FALSE, cm.colors(NumberOfColors))}
		OrderedPoints = order.single(1-DissimilarityMatrix)
		OrderedColors = colors1[OrderedPoints, OrderedPoints]
		PlotParam <- par(mfrow=c(1,2), pty="s")
		MatrixLabels <- colnames(DissimilarityMatrix)
		plotcolors(colors1, main=Title[1])
		plotcolors(OrderedColors, main=Title[2])
		par(PlotParam)
	}

	# Plot the dissimilarity matrices
	for(k in 2:maxColors){
		ColorDissimilarity(DissimilarityMatrix, NumberOfColors=k, byrank, Title)
		readkey()
	}
}

################################################################################
# SeriationPlots: Plot the seriationized dissimilarity matrix from k=2 to maxK #
################################################################################
require(seriation)
SeriationPlots <- function(DissimilarityMatrix, PlotTitle, maxK){
	for(k in 2:maxK){
		ClusterLabels <- kmeans(DissimilarityMatrix, k)$cluster
		dissplot(as.dist(DissimilarityMatrix), labels=ClusterLabels, options = list(key="TRUE", cluster_labels="TRUE", lines=TRUE, silhouettes=TRUE, main=PlotTitle))
		readkey()
	}
}

########################################
# Heatmap Plot of Dissimilarity Matrix #
########################################
HeatmapPlot <- function(DissimilarityMatrix) {
	require(plotly)
	plot_ly(z = DissimilarityMatrix, colorscale = "Hot", type = "heatmap")
}

############################################################################
# 2d scatterplots using Multidimensional Scaling on dissimilarity matrices #
############################################################################
ScatterMDS2D <- function(DissimilarityMatrix){
	require(MASS)
	fit <- cmdscale(as.dist(DissimilarityMatrix), k=2) # k is the number of dim
	x <- fit[,1]
	y <- fit[,2]
	plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="2D Metric MDS", type="p")
	text(x, y, labels = row.names(DissimilarityMatrix), cex=.5) 
}

###################################################################################################
# Interactive 3d scatterplot using Multidimensional Scaling on dissimilarity matrices with plotly #
###################################################################################################
ScatterMDS3D <- function(DissimilarityMatrix){
	library(MASS)
	require(plotly)
	fit <- cmdscale(as.dist(DissimilarityMatrix), k=3) # k is the number of dim
	x <- scale(fit[,1])
	y <- scale(fit[,2])
	z <- scale(fit[,3])
	df <- setNames(data.frame(scale(fit)), c("x", "y", "z"))
	plot_ly(df, x = x, y = y, z = z, type = "scatter3d", mode="markers", marker=list("size"=2,"opacity"=0.75))
}

####################################################################################################
# Interactive 3d surface plot using Multidimensional Scaling on dissimilarity matrices with plotly #
####################################################################################################
SurfaceMDS3D <- function(DissimilarityMatrix){
	library(MASS)
	require(plotly)
	fit <- cmdscale(as.dist(DissimilarityMatrix), k=3) # k is the number of dim
	plot_ly(z = fit, type = "surface")
}

#######################################################################################
# 3d Plot using Multidimensional Scaling on dissimilarity matrices with scatterplot3d #
#######################################################################################
Scatter2MDS3D <- function(DissimilarityMatrix){
	library(MASS)
	require(scatterplot3d)
	fit <- cmdscale(as.dist(DissimilarityMatrix), k=3) # k is the number of dim
	x <- fit[,1]
	y <- fit[,2]
	z <- fit[,3]
	scatterplot3d(x, y, z, highlight.3d = TRUE, col.axis = "blue", col.grid = "lightblue", main = "Helix", pch = 20)
}	


