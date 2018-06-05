#==================================================================
# Load dependencies and call constants
#==================================================================

library(rgdal)
library(raster)
library(rgeos)
library(sp)
options(stringsAsFactors = FALSE)

setwd('..')

imgnum <- 3

img.dest <- paste0('scene', imgnum, 'data', sep="") # for putting data into separate folders for each landsat scene

load(paste0('scene',imgnum, sep=""))

#==================================================================
# Reading imagery
#==================================================================

# # Read Landsat 8 OLI images
# 
# files <- c('imgs\\Bulk Order 910930\\201609_Flood_US', # scene1
#            'imgs\\Bulk Order 912622\\201309_Floods_CO\\LS080320322013090300000000MS00_GO006005004',
#            'imgs\\Bulk Order 912622\\201309_Floods_CO\\LS080310322013082700000000MS00_GO006005004',
#            'imgs\\Bulk Order 912622\\201512_Flood_Midwest_US\\LS080210332016010200000000MS00_GO006005004',
#            'imgs\\Bulk Order 912622\\201311_Floods_TX\\LS080270392013110300000000MS00_GO006005004')
# 
# ## actually only 1-5 because i gave up on reading 6-9 (don't want thermal or pan, 
# # but the files are stored differently from those downloaded from disaster events page)
# # c('imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002007_20170301_20171020_C01_V01_SR',
# # 'imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002008_20170301_20171020_C01_V01_SR',
# # 'imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002008_20180225_20180309_C01_V01_SR',
# # 'imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LE07_CU_002008_20170214_20170928_C01_V01_SR')
# 
# img.names <- list()
# img.meta <- list()
# 
# for(i in 1:length(files)){
#   scene_files <- list.files(files[i],
#                             pattern = '.TIF$|.tif',
#                             full.names = TRUE)[c(1:7,9)] # only gets bands 1-7, not pan or thermal
#   img.names[[i]] <- scene_files
#   img.meta[[i]] <- list.files(files[i], 
#                               pattern = '.txt$',
#                               full.names=TRUE)
# }
# 
# # img.meta <- lapply(img.meta, read.csv) # can't read XML files
# 
# # These are supposed to read XML files but get error: 'attempt to set an attribute on NULL'
# # require(RStoolbox)
# # readMeta('imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002007_20170301_20171020_C01_V01_SR\\LC08_CU_002007_20170301_20171020_C01_V01.xml')
# 
# # stack files
# imgs <- list()
# for(i in 1:length(img.names)){ 
#   scene_stack <- stack(img.names[[i]])
#   scene <- brick(scene_stack)
#   meta <- read.csv(img.meta[[i]])
#   save(scene, meta, file=paste0('scene',i, sep=""))
#   print(paste0('Finished scene', i))
# }
# 
# # devtools::install_github("16EAGLE/getSpatialData")
# # require('getSpatialData')
# # can only get L1 products from Landsat-8
# 
# # plotRGB(scene, r=4,g=3,b=2, stretch='hist')


#==================================================================
# Some constants
#==================================================================

#=================== Lat/long of scene centroids
# For later use in downloading aux data
centroids <- list()
centroids[[1]] <- c(41.746655, -91.41925) # scene1 lat, long
centroids[[2]] <- c(40.319855,-102.72957) # scene2
centroids[[3]] <- c(40.31845, -101.187445) # scene3
centroids[[4]] <- c(38.892855, -86.15885) # scene4
centroids[[5]] <- c(30.295545, -97.865065) # scene5

#=================== urls for JRC permanent water extent tiles
extenturls <- list()
extenturls[[1]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_100W_50N.tif'
extenturls[[2]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_110W_50N.tif'
extenturls[[3]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_110W_50N.tif'
extenturls[[4]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_90W_40N.tif'
extenturls[[5]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_100W_40N.tif'


#==================================================================
# Auxiliary data
#==================================================================

setwd(img.dest) # setwd to the dest file for this scene

#=================== Permanent Water Extent

# use this to find the correct JRC GSW tile 
# https://global-surface-water.appspot.com/download

# Get JRC water surface extent
extenturl <- extenturls[[imgnum]]
download.file(extenturl, 'jrc.tif', mode='wb')

wext <- raster('jrc.tif')
plot(wext)

rasterOptions(chunksize = 1e+04, maxmemory = 1e+06)

# ## Crop and reproject wext to scene
# ext <- extent(scene) # get extent of scene
# ext <- as(ext, "SpatialPolygons") # convert to sp
# ext@proj4string@projargs <- crs(scene)@projargs # set crs of scene extent polygon
# ext <- spTransform(ext, crs(wext)) # convert from scene CRS to wext CRS
# wext.crop <- crop(wext, ext) # crop wext to the scene extent
# writeRaster(wext.crop, "wext_crop_scene.tif", format = "GTiff") # export to local file so we can use with gdalwarp
# 
# # reproject cropped wext to scene CRS
# require(gdalUtils)
# gdalwarp('wext_crop_scene.tif', 'wext_warp_scene.tif', s_srs = crs(wext)@projargs, t_srs = crs(scene)@projargs)
# 
# wext.proj <- raster('wext_warp_scene.tif') # need to remove the extra .tif

# # trying reprojecting and cropping with just gdalwarp
# ext <- extent(scene) # get extent of scene
# ext <- as(ext, "SpatialPolygons") # convert to sp
# ext@proj4string@projargs <- crs(scene)@projargs # set crs of scene extent polygon
# ext <- spTransform(ext, crs(wext)) # convert from scene CRS to wext CRS
# ext <- as(ext, "SpatialPolygonsDataFrame") # for some reason have to do this again for writeOGR
# writeOGR(ext, getwd(), 'ext', drive='ESRI Shapefile') # save extent shapefile
# writeRaster(wext.crop, "wext_crop_scene.tif", format = "GTiff") # export to local file so we can use with gdalwarp
# gdalwarp('jrc.tif', 'wext_warp_scene.tif', # reproject AND crop wext to extent
#          cutline = 'ext.shp',
#          crop_to_cutline = TRUE,
#          s_srs = crs(wext)@projargs, 
#          t_srs = crs(scene)@projargs)
# 
# wext.proj <- raster('wext_warp_scene.tif')

# what about just using align raster?

input <- paste0('..\\', img.names[[imgnum]][1], sep="") # get a raster from the scene
align_rasters('jrc.tif',input, dstfile = 'wext_align.tif')
wext.proj <- raster('wext_align.tif')

#==================================================================
# Auxiliary Data
#==================================================================

#===================  Distance from water bodies
wext.crop <- crop(wext.proj, scene)
wext.crop[wext.crop==0] <- NA # set 0 values to NA
# wext.dist <- distance(wext.crop)

#=================== Get DEM

# require(elevatr)
## get bounding box
# can't figure out zoom level with this
# dem <- get_aws_terrain(bbox(mndwi.clouds), z=5, crs(mndwi.clouds)@projargs)
# plot(dem)
# plot(z, add=TRUE)

# This downloads the srtm that covers the given lat/lon, but more tiles may be necessary to cover the full landsat scene
srtm <- getData('SRTM', lat=centroids[[imgnum]][1], lon=centroids[[imgnum]][2]) # centroid of scene
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


#=================== Distance from Rivers (Area)
require(FedData)
nhd <- get_nhd(template = ext, label='scene')
# nhd.rast <- raster(nrow=dim(scene)[1], ncol=dim(scene)[2])
# extent(nhd.rast) <- extent(nhd$Area)
# nhd.rast <- rasterize(nhd$Area, nhd.rast)


blank <- matrix(nrow=dim(scene[[1]])[1], ncol=dim(scene[[1]])[2])# create empty matrix to burn vector into
blank <- raster(blank)
crs(blank) <- crs(nhd$Area)
extent(blank) <- extent(nhd$Area)
writeRaster(blank, 'blank.tif', format="GTiff", overwrite=TRUE)

nhd.input <- 'EXTRACTIONS\\scene\\NHD\\scene_NHD_Area.shp'
gdal_rasterize(src_datasource = nhd.input, 
               dst_filename = 'blank.tif',
               burn=1)

input <- paste0('..\\', img.names[[imgnum]][1], sep="") # get a raster from the scene
align_rasters('blank.tif',input, dstfile = 'nhd_align.tif')

nhd.rast <- raster('nhd_align.tif')

nhd.dist <- raster('nhd_dist_arc.tif') # compute distance in ArcMap


# reproject raster
# require(gdalUtils)
# writeRaster(nhd.rast, "nhd_rast.tif", format = "GTiff") # export to local file so we can use with gdalwarp
# gdalwarp('nhd_rast.tif', 'nhd_rast_warp.tif',  s_srs = crs(nhd.rast)@projargs, t_srs = crs(scene)@projargs) # reproject 
# nhd.rast.proj <- raster('nhd_rast_warp.tif') # dimensions are slightly different from scene, need to trim

nhd.rast.proj <- crop(nhd.rast.proj, scene)
nhd.dist <- distance(nhd.rast.proj)
 # crop because dims changed after gdalwarp

# try using gdal_grid for distance
nhd.input <- 'EXTRACTIONS\\scene\\NHD\\scene_NHD_Area.shp'
gdal_grid(src_datasource = nhd.input,
          txe = c(extent(scene)[1], extent(scene)[2]),
          tye = c(extent(scene)[3], extent(scene)[4]),
          dst_filename = 'nhd_area_dist.tif',
          a = 'average_distance'
          )
test <- raster('nhd_area_dist.tif')

#=================== Distance from Rivers (Flowlines)
nhd.rast.f <- raster(nrow=dim(scene)[1], ncol=dim(scene)[2])
extent(nhd.rast.f) <- extent(nhd$Flowline)
nhd.rast.f <- rasterize(nhd$Flowline, nhd.rast.f)

# could also use gdal_rasterize?

# reproject raster
require(gdalUtils)
writeRaster(nhd.rast.f, "nhd_rast_f.tif", format = "GTiff") # export to local file so we can use with gdalwarp
gdalwarp('nhd_rast_f.tif', 'nhd_rast_f_warp.tif',  s_srs = crs(nhd.rast)@projargs, t_srs = crs(scene)@projargs)# reproject
nhd.rast.f.proj <- raster('nhd_rast_f_warp.tif') # dimensions are slightly different from scene, need to trim

nhd.rast.f.proj <- crop(nhd.rast.f.proj, scene) # crop because dims changed after gdalwarp
nhd.dist.f <- distance(nhd.rast.f.proj)



#==================================================================
# Detect floodwaters
#==================================================================

#=================== Calculate MNDWI 
# Modified NDWI
# (green - SWIR2) / (green + SWIR2)
mndwi <- (scene[[3]] - scene[[7]]) / (scene[[3]] + scene[[7]])

# Use Otsu histogram thresholding to find cutoff in MNDWI between flood and no flood (maximize variance)
# source("https://bioconductor.org/biocLite.R")
# biocLite("EBImage")
require(EBImage)

# find water vs. non-water threshold using otsu on full image
mndwi.mat <- matrix(getValues(mndwi), nrow=dim(mndwi)[1], ncol=dim(mndwi)[2], byrow=TRUE)
mndwi.otsu <- otsu(mndwi.mat, range=c(-1,1))
water <- mndwi
water[water > mndwi.otsu] <- 1
water[water <= mndwi.otsu] <- 0

water <- reclassify(water, c(mndwi.otsu, max(getValues(water), na.rm=T), 1))
water <- reclassify(water, c(min(getValues(water),na.rm=T), mndwi.otsu, 0))


#==================================================================
# Processing for classifier
#==================================================================
# At this point, all the data will be put into one huge matrix for training in Python

xmin = extent(scene)[1]
xmax = extent(scene)[2]
ymin = extent(scene)[3]
ymax = extent(scene)[4]

# get random samples within image for training and testing
seed=1
n.samples <- 5000
x <- c(runif(n.samples, xmin, xmax))
y <- c(runif(n.samples, ymin, ymax))
samples <- cbind(x,y)

# extract values from each raster
# ... is all input rasters, pts is matrix of sampling point coordinates
test <- function(..., pts){
  rasts <- list(...)
  vals <- lapply(rasts, function(x) extract(x, pts))
  do.call(cbind, vals)
}

rast.mat <- test(slope, aspect, nhd.dist, nlcd.proj, water, pts=samples)
rast.mat <- cbind(samples,rast.mat) # add XY as well

rast.mat <- rast.mat[complete.cases(rast.mat),] # remove rows with NA
write.csv(rast.mat, "rast_mat.csv", row.names=FALSE)


#==================================================================
# Generate clouds
#==================================================================
# Generates structured noise using gaussian random fields

ext <- extent(scene)
ext <- as(ext, "SpatialPolygons") # convert to sp
p<-ext@polygons[[1]]@Polygons[[1]]@coords # get coords of extent box
simdim <- min(dim(scene)[1:2]) # get min dimensions of mndwi, either nrow or ncol, which we will square to make a regular matrix

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


#==================================================================
# Processing for cloudy pixel recovery
#==================================================================
# We want to see how well we can actually guess what is under clouds, so let's create a separate test set of clouded out pixels

# Masking out optical data
mndwi.clouds <- mask(mndwi, z) # mask out clouds in mndwi
scene.clouds <- mask(scene, z) # mask out clouds in multiband scene

#==================================================================
# Save and delete
#==================================================================
# These are big files, so save the ones we might want later as .rda files and delete the rest

#=================== Save as RDA
save(wext.proj, wext.dist, file='wext')
save(srtm.proj, slope, aspect, file='srtm_data')
save(nlcd.proj, file='nlcd')
save(nhd.rast.proj, nhd.dist, nhd.rast.f.proj, nhd.dist.f, file='nhd')
save(mndwi, water, file='water')
save(ppad, z, mndwi.clouds, scene.clouds, file='clouds')

#=================== Delete the rest

# let's try zipping them first, then we can just put them on the U-Drive

bye.files <- c('jrc.tif', 
         'wext_crop_scene.tif', 
         'wext_warp_scene.tif', 
         'srtm_crop_scene.tif', 
         'srtm_warp_scene.tif', 
         'slope.tif', 
         'aspect.tif', 
         'nlcd_scene.tif', 
         'nlcd_warp_scene.tif', 
         'nhd_rast.tif', 
         'nhd_rast_warp.tif', 
         'nhd_rast_f.tif', 
         'nhd_rast_f_warp.tif')

# bye.files <- list.files(getwd())


# zip(zipfile = paste0('scene',imgnum, 'aux_data', sep=""), files=bye.files)
# bye.folders <- c('RAW', 'EXTRACTIONS')

# file.remove(bye.files)
