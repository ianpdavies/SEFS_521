#==================================================================
# Creating plots for final presentation and paper




#==================================================================
# Load dependencies and call constants
#==================================================================

require(ggplot2)
require(raster)
require(rasterVis)

setwd('..')
imgnum <- 1
img.dest <- paste0('scene', imgnum, 'data', sep="") # for putting data into separate folders for each landsat scene
load(paste0('scene',imgnum, sep=""))
setwd(img.dest) # setwd to the dest file for this scene

load('wext')
load('srtm_data')
load('nlcd')
load('nhd')
load('water')
load('cloud_masks')

#==================================================================
# NDWI and MNDWI for plotting in ENVI
#==================================================================

ndwi <- (scene[[3]] - scene[[4]]) / (scene[[3]] + scene[[4]])
writeRaster(ndwi, 'ndwi.tif', format='GTiff')

mndwi <- (scene[[3]] - scene[[7]]) / (scene[[3]] + scene[[7]])
writeRaster(mndwi, 'mndwi.tif', format='GTiff')

#==================================================================
# Study areas
#==================================================================

load('..\\scene1')
plotRGB(scene, r=4, g=3, b=2, stretch='hist', colNA='white', axes=FALSE)

#==================================================================
# Otsu histogram
#==================================================================

require(EBImage)
mndwi.mat <- matrix(getValues(mndwi), nrow=dim(mndwi)[1], ncol=dim(mndwi)[2], byrow=TRUE)
mndwi.otsu <- otsu(mndwi.mat, range=c(-1,1))
hist(mndwi.mat, xlab='MNDWI Value', ylab='Frequency', col='grey', xlim=c(-0.5,0.5), breaks=200, main='')
abline(v=mndwi.otsu,col="red", lwd=2)

#==================================================================
# Flood and non-flood after Otsu thresholding
#==================================================================

plot(water, legend=FALSE,  axes=FALSE, box=FALSE)

#==================================================================
# Cloud simulation and masking
#==================================================================

water.clouds <- mask(water, z)
water.noclouds <- mask(water, y)

breakpoints <- c(-1, 0, 1, NA)
colors <- c('pink','blue', 'blue','blue')

plot(water.clouds, legend=FALSE, box=FALSE, axes=FALSE, colNA='black', breaks=breakpoints, col=colors)
plot(water.noclouds, legend=FALSE, box=FALSE, axes=FALSE, colNA='black', breaks=breakpoints, col=colors)

#==================================================================
# Auxiliary variables
#==================================================================

require(RColorBrewer)

colors <- brewer.pal(n = 10, name = "OrRd")
plot(nhd.dist, axes=FALSE, box=FALSE,col=colors)

plot(nlcd.proj, main='Land cover',  axes=FALSE, box=FALSE)

colors <- brewer.pal(n = 10, name = "RdYlGn")
plot(slope, col=colors, axes=FALSE, box=FALSE)

colors <- brewer.pal(n = 5, name = "Blues")
plot(aspect, col=colors, axes=FALSE, box=FALSE)


#==================================================================
# Cloud prediction
#==================================================================

preds <- read.csv('cloud_preds.csv')
clouds <- read.csv('rast_mat_clouds.csv')
clouds <- clouds[complete.cases(clouds),]

preds_acc <- cbind(clouds, preds)
colnames(preds_acc) <- c('x','y','slope','aspect','dist','lulc','actual', 'delete','pred')
preds_acc$acc <- ifelse(preds_acc$actual==preds_acc$pred, 'green','red')

dev.off()
plot(water.noclouds)
points(x=preds_acc[,1], y=preds_acc[,2], col=preds_acc$acc, pch=20, add=TRUE)

uncloud.samp <- read.csv('uncloud_samp.csv')
plot(x=uncloud.samp[,1], y=uncloud.samp[,2], pch=20)


require(gridExtra)
layout(matrix(1:6,ncol=2))
image(mask.rasts[[1]])
image(mask.rasts[[2]])
image(mask.rasts[[3]])
image(mask.rasts[[4]])
image(mask.rasts[[5]])
image(mask.rasts[[1]])
