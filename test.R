require(rgdal)
require(raster)
require(rgeos)
require(sp)
require(gdalUtils)
require(FedData)
require(data.table) # for `fwrite` to export matrices
options(stringsAsFactors = FALSE)

setwd('..')
imgnum <- 1 # This determines which scene we are looking at
img.dest <- paste0('scene', imgnum, 'data', sep="") # for putting data into separate folders for each landsat scene
load(paste0('scene',imgnum, sep=""))

# # Read Landsat 8 OLI images
files <- c('imgs\\Bulk Order 910930\\201609_Flood_US', # scene1
           'imgs\\Bulk Order 912622\\201309_Floods_CO\\LS080320322013090300000000MS00_GO006005004',
           'imgs\\Bulk Order 912622\\201309_Floods_CO\\LS080310322013082700000000MS00_GO006005004',
           'imgs\\Bulk Order 912622\\201512_Flood_Midwest_US\\LS080210332016010200000000MS00_GO006005004',
           'imgs\\Bulk Order 912622\\201311_Floods_TX\\LS080270392013110300000000MS00_GO006005004')

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

input <- paste0('..\\', img.names[[imgnum]][1], sep="") # the file location of the scene raster, used for GDAL operations

setwd(img.dest)
load('wext', verbose=TRUE)
load('srtm_data', verbose=TRUE)
load('nlcd', verbose=TRUE)
load('nhd', verbose=TRUE)
load('water', verbose=TRUE)
load('cloud_masks', verbose=TRUE)

#==================================================================
# Extracting values from cloudy and non-cloudy pixels
#==================================================================

# Need to increase memory usage for masking. Benchmark gives me a little over 1 minute
rasterOptions(chunksize = 1e+08, maxmemory = 1e+09)

masker <- function(..., m){ # ... is all input rasters, z is raster of
  rasts <- list(...)
  mask.rasts <- lapply(rasts, function(x) mask(x, m)) # will differently sized rasters be a problem?
  mask.vals <- lapply(mask.rasts, values)
  do.call(cbind, mask.vals)
}

# Create a raster with cellnumbers so we can reassemble cells into an image of predictions
cellnums <- slope
cellnums <- setValues(cellnums, values = seq(1, length(cellnums), 1))
writeRaster(cellnums, 'cellnums.tif', format='GTiff', overwrite=TRUE)
cellnums <- raster('cellnums.tif')

#=================== Clouded pixels
rast.mat.clouds <- masker(cellnums, slope, aspect, nhd.dist, nlcd.proj, water, m=z)
rast.mat.clouds <- cbind(xyFromCell(z, 1:length(z)), # add XY as well from cell numbers of z
                         rast.mat.clouds) # add data
rast.mat.clouds <- rast.mat.clouds[complete.cases(rast.mat.clouds),] # remove rows with NA
colnames(rast.mat.clouds) <- c('x','y','cellnums','slope','aspect','nhd_dist', 'nlcd', 'water') # give proper col names
fwrite(data.frame(rast.mat.clouds), "rast_mat_clouds.csv", row.names=FALSE) # save matrix as csv for use in python

#=================== Unclouded pixels
rast.mat.unclouds <- masker(cellnums, slope, aspect, nhd.dist, nlcd.proj, water, m=y)
rast.mat.unclouds <- cbind(xyFromCell(y, 1:length(y)), 
                           rast.mat.unclouds)
rast.mat.unclouds <- rast.mat.unclouds[complete.cases(rast.mat.unclouds),]
colnames(rast.mat.unclouds) <- c('x','y','cellnums', 'slope','aspect','nhd_dist', 'nlcd', 'water')
fwrite(data.frame(rast.mat.unclouds), "rast_mat_unclouds.csv", row.names=FALSE)

 
