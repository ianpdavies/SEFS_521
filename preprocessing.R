#==================================================================
# Load dependencies and call constants
#==================================================================

require(rgdal)
require(raster)
require(rgeos)
require(sp)
require(gdalUtils)
require(FedData)
require(data.table) # for `fwrite` to export matrices
options(stringsAsFactors = FALSE)

setwd('..')

imgnum <- 1

img.dest <- paste0('scene', imgnum, 'data', sep="") # for putting data into separate folders for each landsat scene

load(paste0('scene',imgnum, sep=""))

#==================================================================
# Reading imagery
#==================================================================

# # Read Landsat 8 OLI images
# 
files <- c('imgs\\Bulk Order 910930\\201609_Flood_US', # scene1
           'imgs\\Bulk Order 912622\\201309_Floods_CO\\LS080320322013090300000000MS00_GO006005004',
           'imgs\\Bulk Order 912622\\201309_Floods_CO\\LS080310322013082700000000MS00_GO006005004',
           'imgs\\Bulk Order 912622\\201512_Flood_Midwest_US\\LS080210332016010200000000MS00_GO006005004',
           'imgs\\Bulk Order 912622\\201311_Floods_TX\\LS080270392013110300000000MS00_GO006005004')

# ## actually only 1-5 because i gave up on reading 6-9 (don't want thermal or pan, 
# # but the files are stored differently from those downloaded from disaster events page)
# # c('imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002007_20170301_20171020_C01_V01_SR',
# # 'imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002008_20170301_20171020_C01_V01_SR',
# # 'imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LC08_CU_002008_20180225_20180309_C01_V01_SR',
# # 'imgs\\Bulk Order 910931\\U.S. Landsat 4-8 ARD\\LE07_CU_002008_20170214_20170928_C01_V01_SR')
# 
img.names <- list()
img.meta <- list()

# gets all layer names for each scene
for(i in 1:length(files)){
  scene_files <- list.files(files[i],
                            pattern = '.TIF$|.tif',
                            full.names = TRUE)[c(1:7,9)] # only gets bands 1-7, not pan or thermal
  img.names[[i]] <- scene_files
  img.meta[[i]] <- list.files(files[i],
                              pattern = '.txt$',
                              full.names=TRUE)
}
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

# urls for JRC permanent water extent tiles
extenturls <- list()
extenturls[[1]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_100W_50N.tif'
extenturls[[2]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_110W_50N.tif'
extenturls[[3]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_110W_50N.tif'
extenturls[[4]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_90W_40N.tif'
extenturls[[5]] <- 'https://storage.googleapis.com/global-surface-water/downloads/extent/extent_100W_40N.tif'

rasterOptions(chunksize = 1e+04, maxmemory = 1e+06)

input <- paste0('..\\', img.names[[imgnum]][1], sep="") # the file location of the scene raster, used for GDAL operations

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
align_rasters('jrc.tif',input, dstfile = 'wext_align.tif') # crop and reproject to scene specifications
wext.proj <- raster('wext_align.tif')

# Create layer of distances from permanent water extent
# ...

#=================== Get DEM

# Find all tiles within scene extent using the SRTM tile grid as a reference
# from https://www.gis-blog.com/download-srtm-for-an-entire-country/
srtm.grid <- shapefile("../srtm/tiles.shp")
scene.shp <- as(extent(scene), 'SpatialPolygons') # create sp from scene extent
crs(scene.shp) <- crs(scene) # add proj
scene.shp <- spTransform(scene.shp, crs(srtm.grid)) # reproject to scene projection
intersects <- gIntersects(scene.shp, srtm.grid, byid=T) # find SRTM tiles that intersect scene
tiles <- srtm.grid[intersects[,1],] # finds tiles that intersect scene, by ID

# Download tiles
srtm_list  <- list()
for(i in 1:length(tiles)) {
  lon <- extent(tiles[i,])[1]  + (extent(tiles[i,])[2] - extent(tiles[i,])[1]) / 2
  lat <- extent(tiles[i,])[3]  + (extent(tiles[i,])[4] - extent(tiles[i,])[3]) / 2
  
  tile <- getData('SRTM', 
                  lon=lon, 
                  lat=lat)
  
  srtm_list[[i]] <- tile
}

# Mosaic tiles
srtm_list$fun <- mean # function for mosaic to compute cell values where layers overlap
srtm_mosaic   <- do.call(mosaic, srtm_list) # mosaic

# Export mosaic
writeRaster(srtm_mosaic, 'srtm_mosaic.tif', format="GTiff", overwrite=TRUE)

#Crop tiles to scene extent
align_rasters('srtm_mosaic.tif', input, dstfile = 'srtm_align.tif')
srtm.proj <- raster('srtm_align.tif')

#=================== calculate slope
gdaldem(mode="slope", input_dem='srtm_align.tif', p=TRUE, output='slope.tif', output_Raster=TRUE, verbose=TRUE)
align_rasters('slope.tif', input, dstfile = 'slope_align.tif', overwrite=TRUE)
slope <- raster('slope_align.tif')

#=================== calculate aspect
gdaldem(mode="aspect", input_dem='srtm_align.tif', p=TRUE, output='aspect.tif', output_Raster=TRUE, verbose=TRUE)
align_rasters('aspect.tif', input, dstfile = 'aspect_align.tif', overwrite=TRUE)
aspect <- raster('aspect_align.tif')

#=================== LULC
require(FedData)
nlcd <- get_nlcd(template = ext, label='scene', year=2011, dataset="landcover")
writeRaster(nlcd, "nlcd_scene.tif", format = "GTiff") # export to local file so we can use with gdalwarp
input <- paste0('..\\', img.names[[imgnum]][1], sep="") # get a raster from the scene
align_rasters('nlcd_scene.tif', input, dstfile = 'nlcd_align.tif')
nlcd.proj <- raster('nlcd_align.tif')

# consider cleaning up the NLCD a bit using raster::focal()

#=================== Distance from Rivers (Area)
nhd <- get_nhd(template = ext, label='scene') # download NHD data (area, flowlines, waterbodies - no way to pick and choose)
nhd <- list(Area = readOGR('EXTRACTIONS\\scene\\NHD\\scene_NHD_Area.shp'),
            Flowline=readOGR('EXTRACTIONS\\scene\\NHD\\scene_NHD_Flowline.shp'),
            Line=readOGR('EXTRACTIONS\\scene\\NHD\\scene_NHD_Line.shp'),
            Waterbody=readOGR('EXTRACTIONS\\scene\\NHD\\scene_NHD_Waterbody.shp')) # just in case nhd is already downloaded

blank <- matrix(nrow=dim(scene[[1]])[1], ncol=dim(scene[[1]])[2])# create empty matrix to burn vector into
blank <- raster(blank)

crs(blank) <- crs(nhd$Area) 
extent(blank) <- extent(nhd$Area)
writeRaster(blank, 'blank.tif', format="GTiff", overwrite=TRUE)

nhd.input <- 'EXTRACTIONS\\scene\\NHD\\scene_NHD_Area.shp'
gdal_rasterize(src_datasource = nhd.input, # burns NHD vector values of 1 into blank raster
               dst_filename = 'blank.tif',
               burn=1)

input <- paste0('..\\', img.names[[imgnum]][1], sep="") # get a raster from the scene
align_rasters('blank.tif',input, dstfile = 'nhd_align.tif') # crop and reproject to scene specifications

nhd.dist <- raster('nhd_dist_arc.tif') # compute distance in ArcMap because distance() takes WAY too long

#=================== Distance from Rivers (Flowlines)
blank <- matrix(nrow=dim(scene[[1]])[1], ncol=dim(scene[[1]])[2])# create empty matrix to burn vector into
blank <- raster(blank) 
crs(blank) <- crs(nhd$Flowline) 
extent(blank) <- extent(nhd$Flowline)
writeRaster(blank, 'blank.tif', format="GTiff", overwrite=TRUE)

nhd.input <- 'EXTRACTIONS\\scene\\NHD\\scene_NHD_Flowline.shp'
gdal_rasterize(src_datasource = nhd.input, # burns NHD vector values of 1 into blank raster
               dst_filename = 'blank.tif',
               burn=1)

input <- paste0('..\\', img.names[[imgnum]][1], sep="") # get a raster from the scene
align_rasters('blank.tif',input, dstfile = 'nhd_f_align.tif') # crop and reproject to scene specifications

nhd.dist.f <- raster('nhd_dist_f_arc.tif') # compute distance in ArcMap EuclideanDistance because distance() takes WAY too long

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
# water[water > mndwi.otsu] <- 1
# water[water <= mndwi.otsu] <- 0

water <- reclassify(water, c(mndwi.otsu, max(getValues(water), na.rm=T), 1))
water <- reclassify(water, c(min(getValues(water),na.rm=T), mndwi.otsu, 0))
# writeRaster(water, 'water.tif')
# 
# input <- paste0('..\\', img.names[[imgnum]][1], sep="") # get a raster from the scene
# align_rasters('water.tif',input, dstfile = 'water.tif') # crop and reproject to scene specifications
# water <- raster('water.tif')

#==================================================================
# Processing for classifier
#==================================================================
# At this point, all the data will be put into one huge matrix for training in Python

# xmin = extent(scene)[1]
# xmax = extent(scene)[2]
# ymin = extent(scene)[3]
# ymax = extent(scene)[4]
# 
# # get random samples within image for training and testing
# seed=1
# n.samples <- 10000
# x <- c(runif(n.samples, xmin, xmax))
# y <- c(runif(n.samples, ymin, ymax))
# samples <- cbind(x,y)
# 
# # extract values from each raster
# # ... is all input rasters, pts is matrix of sampling point coordinates
# test <- function(..., pts){
#   rasts <- list(...)
#   vals <- lapply(rasts, function(x) extract(x, pts))
#   do.call(cbind, vals)
# }
# 
# rast.mat <- test(slope, aspect, nhd.dist, nlcd.proj, water, pts=samples)
# rast.mat <- cbind(samples,rast.mat) # add XY as well
# 
# rast.mat <- rast.mat[complete.cases(rast.mat),] # remove rows with NA
# write.csv(rast.mat, "rast_mat.csv", row.names=FALSE)

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

# Now create rasters from the grf output matrix - one for clouded pixels, one for unclouded pixels (the inverse)
ppad_cloud <- ifelse(ppad>0.3, 1, NA) # reclassify matrix with arbitrary cutoff (< cutoff is unclouded)
z <- raster(nrows=nrow(ppad_cloud), ncols=ncol(ppad_cloud), crs=crs(scene))
extent(z) <- extent(scene)
z <- setValues(z, values=ppad_cloud)

writeRaster(z, 'cloud_vals.tif') # save raster and read it, so the values are stored in file instead of in-memory (faster)
z <- raster('cloud_vals.tif')

ppad_uncloud <- ifelse(ppad<0.3, 1, NA) # reclassify matrix with arbitrary cutoff (< cutoff is clouded)
y <- raster(nrows=nrow(ppad_uncloud), ncols=ncol(ppad_uncloud), crs=crs(scene))
extent(y) <- extent(scene)
y <- setValues(y, values=ppad_uncloud)

writeRaster(z, 'uncloud_vals.tif')
y <- raster('uncloud_vals.tif')

rm(ppad, ppad_cloud, ppad_uncloud)
#==================================================================
# Extracting values from cloudy and non-cloudy pixels
#==================================================================
# create separate masked rasters for unclouded pixels (train) and cloudy pixels (test)

# extract clouded or unclouded pixel values from each raster
masker <- function(..., mask){ # ... is all input rasters, z is raster of
  start.time <- Sys.time()
  rasts <- list(...)
  mask.rasts <- lapply(rasts, function(x) mask(x, z)) # will differently sized rasters be a problem?
  mask.vals <- lapply(mask.rasts, values)
  do.call(cbind, mask.vals)
  paste0("Runtime is", Sys.time() - start.time)
}

#=================== Clouded pixels
rast.mat.clouds <- masker(slope, aspect, nhd.dist, nlcd.proj, water, mask=z)
rast.mat.clouds <- cbind(xyFromCell(z, 1:length(z)), rast.mat.clouds) # add XY as well from cell numbers of z
rast.mat.clouds <- rast.mat.clouds[complete.cases(rast.mat.clouds),] # remove rows with NA
fwrite(data.frame(rast.mat.clouds), "rast_mat_clouds.csv", row.names=FALSE)


#=================== Unclouded pixels
rast.mat.unclouds <- masker(slope, aspect, nhd.dist, nlcd.proj, water, mask=y)
rast.mat.unclouds <- cbind(xyFromCell(y, 1:length(y)), rast.mat.unclouds) # add XY as well from cell numbers of z
rast.mat.unclouds <- rast.mat.unclouds[complete.cases(rast.mat.unclouds),] # remove rows with NA
fwrite(data.frame(rast.mat.unclouds), "rast_mat_unclouds.csv", row.names=FALSE)

#==================================================================
# Save and delete
#==================================================================
# These are big files, so save the ones we might want later as .rda files and delete the rest

#=================== Save as RDA
# save(wext.proj, file='wext')
# save(srtm.proj, slope, aspect, file='srtm_data')
# save(nlcd.proj, file='nlcd')
# save(nhd.rast.proj, nhd.dist, nhd.rast.f.proj, nhd.dist.f, file='nhd')
# save(mndwi, water, file='water')
# save(ppad, z, y, file='cloud_masks')

load('wext', verbose=TRUE)
load('srtm_data', verbose=TRUE)
load('nlcd', verbose=TRUE)
load('nhd', verbose=TRUE)
load('water', verbose=TRUE)
load('cloud_masks', verbose=TRUE)

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
