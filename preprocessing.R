#==================================================================
# Load dependencies
#==================================================================

library(rgdal)
library(raster)
library(rgeos)
library(sp)
options(stringsAsFactors = FALSE)

setwd('C:\\Users\\ipdavies\\Documents\\Floods')

#==================================================================
# Reading imagery
#==================================================================

# Read Landsat 8 OLI images
scene1_files <- list.files('imgs\\Bulk Order 910930\\201609_Flood_US',
                           pattern = '.TIF$',
                           full.names = TRUE)[c(1:7,9)]

scene1_meta <- list.files('imgs\\Bulk Order 910930\\201609_Flood_US',
                          pattern = '.txt$',
                          full.names=TRUE)

scene1_meta <- read.csv(scene1_meta)

# stack files
scene1_stack <- stack(scene1_files)
scene1 <- brick(scene1_stack)

#==================================================================
# Crop and reproject data
#==================================================================

# use this to find the correct JRC GSW tile 
# https://global-surface-water.appspot.com/download

# Get JRC water surface extent
extenturl <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_100W_50N.tif'
download.file(extenturl, 'jrc.tif', mode='wb')
wext <- raster('jrc.tif')
plot(wext)

rasterOptions(chunksize = 1e+04, maxmemory = 1e+06)

## Crop and reproject wext to scene
ext <- extent(scene1) #get extent of scene
ext <- as(ext, "SpatialPolygons") # convert to sp
ext@proj4string@projargs <- crs(scene1)@projargs # set crs of scene extent polygon
ext <- spTransform(ext, crs(wext)) # convert from scene CRS to wext CRS
wext.crop <- crop(wext, ext) # crop wext to the scene extent
writeRaster(wext.crop, "wext_crop_scene1.tif", format = "GTiff") # export to local file so we can use with gdalwarp

# reproject cropped wext to scene CRS
require(gdalUtils)
gdalwarp('wext_crop_scene1.tif', 'wext_warp_scene1.tif.tif', s_srs = crs(wext)@projargs, t_srs = crs(scene1)@projargs)

wext.proj <- raster('wext_warp_scene1.tif.tif') # need to remove the extra .tif


# Maybe use this to find more landsat data?
# devtools::install_github("16EAGLE/getSpatialData")

# maybe use sparse matrices when masking?

#==================================================================
# Detect flood waters
#==================================================================

# calculate modified NDWI
# mndwi <- (scene1[[3]] - scene1[[7]]) / (scene1[[3]] + scene1[[7]])

# save(scene1, scene1_meta, file="scene1")

load("scene1")

plot(mndwi)

#==================================================================
# Generate cloudy pixels
#==================================================================

# generates perlin noise image
# function from: https://stackoverflow.com/questions/15387328/realistic-simulated-elevation-data-in-r-perlin-noise

perlin_noise <- function( 
  n = 5,   m = 7,    # Size of the grid for the vector field
  N = 100, M = 100   # Dimension of the image
) {
  # For each point on this n*m grid, choose a unit 1 vector
  vector_field <- apply(
    array( rnorm( 2 * n * m ), dim = c(2,n,m) ),
    2:3,
    function(u) u / sqrt(sum(u^2))
  )
  f <- function(x,y) {
    # Find the grid cell in which the point (x,y) is
    i <- floor(x)
    j <- floor(y)
    stopifnot( i >= 1 || j >= 1 || i < n || j < m )
    # The 4 vectors, from the vector field, at the vertices of the square
    v1 <- vector_field[,i,j]
    v2 <- vector_field[,i+1,j]
    v3 <- vector_field[,i,j+1]
    v4 <- vector_field[,i+1,j+1]
    # Vectors from the point to the vertices
    u1 <- c(x,y) - c(i,j)
    u2 <- c(x,y) - c(i+1,j)
    u3 <- c(x,y) - c(i,j+1)
    u4 <- c(x,y) - c(i+1,j+1)
    # Scalar products
    a1 <- sum( v1 * u1 )
    a2 <- sum( v2 * u2 )
    a3 <- sum( v3 * u3 )
    a4 <- sum( v4 * u4 )
    # Weighted average of the scalar products
    s <- function(p) 3 * p^2 - 2 * p^3
    p <- s( x - i )
    q <- s( y - j )
    b1 <- (1-p)*a1 + p*a2
    b2 <- (1-p)*a3 + p*a4
    (1-q) * b1 + q * b2
  }
  xs <- seq(from = 1, to = n, length = N+1)[-(N+1)]
  ys <- seq(from = 1, to = m, length = M+1)[-(M+1)]
  outer( xs, ys, Vectorize(f) )
}

# image( perlin_noise() )

a <- .6
k <- 8
m <- perlin_noise(2,2,2^k,2^k)
for( i in 2:k )
  m <- m + a^i * perlin_noise(2^i,2^i,2^k,2^k)
image(m)
m[] <- rank(m) # Histogram equalization
image(m)

# rescale values from 0-255
m <- m/max(m) * 255


##
# alternative method using Gaussian random fields
require(geoR)

sim <- grf(10000, grid="reg", cov.pars=c(0.75, .25))
.Random.seed <- sim$.Random.seed
image(sim, col=gray(seq(1, .1, l=30)))



# Combine clouds with flood scenes using alpha blending
# OR just remove the pixels using a cloud mask for all noise values > X
# https://howtoenmap.jimdo.com/r-remote-sensing/cloud-masking/

#==================================================================
# Auxiliary data
#==================================================================

# Distance from water bodies
# Slope
# Land cover
# include lat/long because flooding will be spatially autocorrelated?

#==================================================================
# Training Classifier
#==================================================================


#==================================================================
# Testing workflow with a smaller image
#==================================================================

# get smaller scene to test
ext <- extent(scene1)/10
scene <- crop(scene1, ext)


#=================== Calculate MNDWI
# (green - SWIR2) / (green + SWIR2)
mndwi <- (scene1[[3]] - scene1[[7]]) / (scene1[[3]] + scene1[[7]])
mndwi <- crop(mndwi, scene)

#### create clouds
ext <- as(ext, "SpatialPolygons") # convert to sp
p<-ext@polygons[[1]]@Polygons[[1]]@coords 

simdim <- min(dim(scene)[1:2]) # get min dimensions of mndwi, either nrow or ncol
require(geoR)
sim <- grf(simdim^2, grid="reg", cov.pars=c(0.75, .25)) # simulate noise
.Random.seed <- sim$.Random.seed # set random seed 
# image(sim, col=gray(seq(1, .1, l=30))) # visualize

# grf only creates regular matrices, so we have to create a blank irregular matrix and fill it with grf output
ppad <- matrix(NA, nrow=dim(scene)[1], ncol=dim(scene)[2]) # create empty matrix of dims of scene
p <- matrix(sim$data, nrow=sqrt(length(sim$data)), byrow=TRUE) # create matrix of cloud values
ppad[1:nrow(p), 1:ncol(p)] <- p # set values of empty matrix to clouds, with NAs padding until cloud extent = scene extent
p <- raster(ppad)
z<-p>0.3 # eliminate some values to have only partially cloudy image

require(gdalUtils)
extent(z) <- extent(scene) # give extent values to clouds
crs(z) <- crs(scene)# set crs of clouds to scene

z[z==0] <- NA # set 0 values to NA for transparency
# plot(mndwi)
# plot(z, add=TRUE)

mndwi.clouds <- mask(mndwi, z) # mask out clouds in mndwi
scene.clouds <- mask(scene, z) # mask out clouds in multiband scene

#===================  Distance from water bodies
wext.crop <- crop(wext.proj, scene)
wext.crop[wext.crop==0] <- NA # set 0 values to NA
wext.dist <- distance(wext.crop)
plot(wext.dist)

#=================== Get DEM

# require(elevatr)
## get bounding box
# can't figure out zoom level with this
# dem <- get_aws_terrain(bbox(mndwi.clouds), z=5, crs(mndwi.clouds)@projargs)
# plot(dem)
# plot(z, add=TRUE)

# This downloads the srtm that covers the given lat/lon, but more tiles may be necessary to cover the full landsat scene
srtm <- getData('SRTM', lon=-91.41925, lat=41.746655) # coords of scene1
ext <- extent(scene) # get extent of scene
ext <- as(ext, "SpatialPolygons") # convert to sp
ext@proj4string@projargs <- crs(scene)@projargs # set crs of scene extent polygon
ext <- spTransform(ext, crs(srtm)) # convert from scene CRS to srtm CRS
srtm.crop <- crop(srtm, ext) # crop srtm to the scene extent
writeRaster(srtm.crop, "srtm_crop_scene.tif", format = "GTiff") # export to local file so we can use with gdalwarp

# reproject cropped srtm to scene CRS
require(gdalUtils)
gdalwarp('srtm_crop_scene.tif', 'srtm_warp_scene.tif', s_srs = crs(srtm)@projargs, t_srs = crs(scene)@projargs)
srtm.proj <- raster('srtm_warp_scene.tif')

#=================== calculate slope
slope <- gdaldem(mode="slope", input_dem='srtm_warp_scene.tif', p=TRUE, output='slope.tif', output_Raster=TRUE, verbose=TRUE)
# need to resample to 30m res, fix dims
# or maybe not if we just  use a matrix instead of a rasterbrick

#=================== calculate aspect
aspect <- gdaldem(mode="aspect", input_dem='srtm_warp_scene.tif', p=TRUE, output='aspect.tif', output_Raster=TRUE, verbose=TRUE)
# need to resample to 30m res, fix dims

#=================== LULC
require(FedData)
nlcd <- get_nlcd(template = ext, label='scene', year=2011, dataset="landcover")
writeRaster(nlcd, "nlcd_scene.tif", format = "GTiff") # export to local file so we can use with gdalwarp
gdalwarp('nlcd_scene.tif', 'nlcd_warp_scene.tif',  s_srs = crs(nlcd)@projargs, t_srs = crs(scene)@projargs) # reproject 
nlcd.proj <- raster('nlcd_warp_scene.tif')

# consider cleaning up the NLCD a bit using raster::focal()

#=================== Load water extent
# code above

wext <- crop(wext.proj, scene)

#=================== Distance from Rivers (Area)
require(FedData)
nhd <- get_nhd(template = ext, label='scene')
nhd.rast <- raster(nrow=dim(scene)[1], ncol=dim(scene)[2])
extent(nhd.rast) <- extent(nhd$`_Area`)
nhd.rast <- rasterize(nhd$`_Area`, nhd.rast)

# reproject raster
require(gdalUtils)
writeRaster(nhd.rast, "nhd_rast.tif", format = "GTiff") # export to local file so we can use with gdalwarp
gdalwarp('nhd_rast.tif', 'nhd_rast_warp.tif',  s_srs = crs(nhd.rast)@projargs, t_srs = crs(scene)@projargs) # reproject 
nhd.rast.proj <- raster('nhd_rast_warp.tif') # dimensions are slightly different from scene, need to trim

nhd.dist <- distance(nhd.rast.proj)
nhd.dist <- crop(nhd.dist, scene) # crop because dims changed after gdalwarp


#=================== Distance from Rivers (Flowlines)
nhd.rast.f <- raster(nrow=dim(scene)[1], ncol=dim(scene)[2])
extent(nhd.rast.f) <- extent(nhd$`_Flowline`)
nhd.rast.f <- rasterize(nhd$`_Flowline`, nhd.rast.f)

# could also use gdal_rasterize?

# reproject raster
require(gdalUtils)
writeRaster(nhd.rast.f, "nhd_rast_f.tif", format = "GTiff") # export to local file so we can use with gdalwarp
gdalwarp('nhd_rast_f.tif', 'nhd_rast_f_warp.tif',  s_srs = crs(nhd.rast)@projargs, t_srs = crs(scene)@projargs)# reproject 
nhd.rast.f.proj <- raster('nhd_rast_f_warp.tif') # dimensions are slightly different from scene, need to trim

nhd.dist.f <- distance(nhd.rast.f.proj)
nhd.dist.f <- crop(nhd.dist.f, scene) # crop because dims changed after gdalwarp

#=================== XY data

# Inundation is spatially autocorrelated. 
xy <- rasterToPoints(scene)[,c(1:2)] # just the coords

#=================== Detect flood waters for entire scene to use as response variable in classifier

# export mndwi raster to ENVI
writeRaster(mndwi, "mndwi.tif", format="GTiff", overwrite=TRUE)

# export scene
writeRaster(scene, 'scene.tif', format='GTiff')

# export ndvi
ndvi <- (scene1[[5]] - scene1[[4]]) / (scene1[[5]] + scene1[[4]])
ndvi <- crop(ndvi, scene)
writeRaster(ndvi, 'ndvi.tif', format='GTiff')
plot(ndvi)

# get training points in ENVI, read them here


# Otsu
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
require(EBImage)

# find water vs. non-water threshold using otsu on full image
mndwi.full <- (scene1[[3]] - scene1[[7]]) / (scene1[[3]] + scene1[[7]])
mndwi.mat <- matrix(getValues(mndwi.full), nrow=dim(mndwi.full)[1], ncol=dim(mndwi.full)[2], byrow=TRUE)
mndwi.otsu <- otsu(mndwi.mat, range=c(-1,1))
water <- mndwi
water[water > mndwi.otsu] <- 1
water[water <= mndwi.otsu] <- 0

#=================== Training classifier

xmin = extent(scene)[1]
xmax = extent(scene)[2]
ymin = extent(scene)[3]
ymax = extent(scene)[4]

# get random samples within image
seed=1
n.samples <- 5000
x <- c(runif(n.samples, xmin, xmax))
y <- c(runif(n.samples, ymin, ymax))
samples <- cbind(x,y)

# extract values from each raster
vals <- extract(mndwi, samples)
vals <- lapply(list(mndwi, slope), function(x) extract(x, samples))
cbind(vals[[1]],vals[[2]])
do.call(cbind,vals)

# ... is all input rasters, pts is matrix of xy for sampling points
test <- function(..., pts){
  rasts <- list(...)
  vals <- lapply(rasts, function(x) extract(x, pts))
  do.call(cbind, vals)
}

rast.mat <- test(slope, aspect, nhd.dist, nlcd.proj, water, pts=samples)
rast.mat <- cbind(samples,rast.mat)

rast.mat <- rast.mat[complete.cases(rast.mat),] # remove rows with NA
write.csv(rast.mat, "rast_mat.csv", row.names=FALSE)

# split into training and testing
n <- nrow(rast.mat)
seed=1
shuffled_mat <- rast.mat[sample(n), ]
train_indices <- 1:round(0.6 * n)
train <- shuffled_mat[train_indices, ]
test_indices <- (round(0.6 * n) + 1):n
test <- shuffled_mat[test_indices, ]

# KNN

# Random Forest
train.pts <- SpatialPointsDataFrame(train[,c(1:2)], train, proj4string=crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs")) # create spatial dataframe with training points

sclass <- superClass(r, trainData=train, responseCol="class",model = "rf", tuneLength = 1) # supervised classification

require(randomForest)


# save(slope,
#      aspect,
#      scene.clouds,
#      scene,
#      p,
#      srtm.proj,
#      wext.proj,
#      ext,
#      mndwi,
#      mndwi.clouds,
#      nlcd.proj,
#      nhd,
#      nhd.dist,
#      p,
#      file="aux_data")

load('aux_data')
