#' calculate_viewmetrics
#' @description The calculate_viewmetrics function is designed to compute a set of
#' configuration metrics based on a given viewshed object and optionally, digital surface
#' models (DSM) and digital terrain models (DTM) for terrain analysis.
#' The function calculates various metrics that describe the visibility characteristics
#' of a landscape from a specific viewpoint.The metrics include:
#' 1. Extent: The total area of the viewshed, calculated as the number of visible grid
#' cells multiplied by the grid resolution
#' 2. Depth: The furthest visible distance within the viewshed from the viewpoint
#' 3. Vdepth: The standard deviation of distances to visible points,
#' providing a measure of the variation in visible distances
#' 4. Horizontal: The total visible horizontal or terrestrial area within the viewshed
#' 5. Relief: The standard deviation of elevations of the visible ground surface
#' 6. Skyline: The standard deviation of the vertical viewscape, including visible
#' canopy and buildings, when specified
#' 7. Number of patches: Visible fragmentation measured by total visible patches
#' with the viewscape
#' 8. Mean shape index: Visible patchiness based on average perimeter-to-area ratio
#' for all viewscape patches (vegetation and building)
#' 9. Edge density: A measure of visible complexity based on the length of
#' patch edges per unit area
#' 10. Patch size: Total average size of a patches over the entire viewscape area
#' 11. Patch density: Visible landscape granularity based on measuring patch density
#' 12. Shannon diversity index: The abundance and evenness of land cover/use in a viewshed
#' 13. Proportion of object: Proportion of a single type of land use or cover in a viewshed
#'
#' @param viewshed Viewshed object.
#' @param dsm Raster, Digital Surface Model for the calculation of
#' @param dtm Raster, Digital Terrain Model
#' @param masks List, a list including rasters of canopy and building footprints.
#' For example of canopy raster, the value for cells without canopy should be 0 and
#' the value for cells with canopy can be any number.
#' @return List
#' @references Tabrizian, P., Baran, P.K., Berkel, D.B., Mitásová, H., & Meentemeyer, R.K. (2020).
#' Modeling restorative potential of urban environments by coupling viewscape analysis of lidar
#' data with experiments in immersive virtual environments. Landscape and Urban Planning, 195, 103704.
#' @import terra
#' @importFrom terra patches
#' @importFrom terra as.polygons
#' @importFrom terra mask
#' @importFrom terra crop
#' @importFrom ForestTools vwf
#' @importFrom ForestTools mcws
#'
#' @examples
#' \donttest{
#' # Load in DSM
#' test_dsm <- terra::rast(system.file("test_dsm.tif",
#'                                     package ="viewscape"))
#' # Load DTM
#' test_dtm <- terra::rast(system.file("test_dtm.tif",
#'                                     package ="viewscape"))
#'
#' # Load canopy raster
#' test_canopy <- terra::rast(system.file("test_canopy.tif",
#'                                        package ="viewscape"))
#'
#' # Load building footprints raster
#' test_building <- terra::rast(system.file("test_building.tif",
#'                                          package ="viewscape"))
#'
#' # Load in the viewpoint
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp",
#'                                           package = "viewscape"))
#'
#' # Compute viewshed
#' output <- viewscape::compute_viewshed(dsm = test_dsm,
#'                                       viewpoints = test_viewpoint,
#'                                       offset_viewpoint = 6, r = 1600,
#'                                       raster=FALSE)
#'
#' # calculate metrics given the viewshed, canopy, and building footprints
#' test_metrics <- calculate_viewmetrics(viewshed=output,
#'                                                  dsm=test_dsm,
#'                                                  dtm=test_dtm,
#'                                                  horizontal_raster=TRUE
#'                                                  )
#' }
#'
#' @export

generate_circle_points <- function(center, distance) {
  # center: numeric vector c(x, y)
  # distance: radius (in same units as center)

  if (length(center) != 2) {
    stop("Center must be a numeric vector of length 2: c(x, y)")
  }

  # Generate 36 angles (every 10 degrees)
  angles <- seq(0, 350, by = 10)

  # Convert to radians
  rad <- angles * pi / 180

  # Compute x/y coordinates
  x <- center[1] + distance * cos(rad)
  y <- center[2] + distance * sin(rad)

  # Return as a data frame of coordinates
  data.frame(id = 1:36, x = x, y = y)
}


calculate_viewmetrics <- function(viewshed, dsm, dtm, masks = list(),
                                  horizontal_raster=FALSE) {
  if (missing(viewshed)) {
    stop("Viewshed object is missing")
  }
  if (missing(dsm) || missing(dtm)) {
    stop("DSM or DTM is missing")
  }
  units <- sf::st_crs(viewshed@crs)$units
  if (units == "ft") {
    error <- 1.6
    minHeight <- 10
  } else if (units %in% c("m", "metre", "meter")){
    error <- 0.5
    minHeight <- 3
  }
  if (isFALSE(terra::crs(dsm, proj = TRUE) == viewshed@crs)) {
    dsm <- terra::project(dsm, y=terra::crs(viewshed@crs))
  }
  if (isFALSE(terra::crs(dtm, proj = TRUE) == viewshed@crs)) {
    dsm <- terra::project(dtm, y=terra::crs(viewshed@crs))
  }
  output <- list()
  ## visiblepoints: tbale with x y cols
  visiblepoints <- filter_invisible(viewshed, FALSE)
  # viewshed raster, ie spatVector with point geometry
  m <- terra::vect(sp::SpatialPoints(visiblepoints))
  terra::crs(m) <- viewshed@crs
  # Mask raster equals visible viewshed
  mask_ <- terra::mask(filter_invisible(viewshed, TRUE), m)
  # get subdsm/dtm
  subdsm <- terra::crop(dsm, terra::ext(mask_))
  subdtm <- terra::crop(dtm, terra::ext(mask_))
  submodel <- subdsm - subdtm

  # Try to extract crowns.
  results <- tryCatch({
    ttops <- ForestTools::vwf(CHM = submodel,
                              winFun = function(x){x * 0.05 + 0.6},
                              minHeight = minHeight)
    crowns <- ForestTools::mcws(treetops = ttops,
                                CHM = submodel,
                                minHeight = minHeight)
    crowns <- terra::patches(crowns, directions=4)
    crowns <- terra::as.polygons(crowns)
    # viewshed patch parameters
    patch_paras <- patch_p(mask_, crowns)
    list(patch_paras=patch_paras,
         x=patch_paras[[6]][,1],
         y=patch_paras[[6]][,2])
    }, error=function(e){
      print(e)
      print('Defaulting to 0 crowns')
      outer_points <- generate_circle_points(c(viewshed@viewpoint[1],
                                  viewshed@viewpoint[2]),
                                distance=max(viewshed@extent[3]-viewshed@extent[1],
                                             viewshed@extent[4]-viewshed@extent[2])/2)
    list(patch_paras=list(0, 0, 0, 0, 0, 0),
         x=outer_points$x,
         y=outer_points$y)
    }
  )
  patch_paras <- results$patch_paras
  x <- results$x
  y <- results$y

  pointnumber <- length(visiblepoints[,1])
  resolution <- viewshed@resolution[1]
  # Number of patches
  # Mean shape index
  # Edge density
  # Patch size
  for (i in 1:5) {
    output[[length(output)+1]] <- patch_paras[[i]]
  }
  names(output) <- c("Nump", "MSI", "ED", "PS", "PD")
  # extent - Total area of the viewshed
  extent <- pointnumber * resolution^2
  output[[length(output)+1]] <- extent
  # depth - Furthest visible distance given the viewscape
  depths <- get_depths(viewshed@viewpoint[1],  # Viewpoint X-coordinate
                       viewshed@viewpoint[2],  # Viewpoint y-coordinate
                       x,  # x-coordinates of
                       y,
                       length(x))
  depths <- depths[!is.na(depths)]
  # depth
  output[[length(output)+1]] <- max(depths)
  # vdepth
  output[[length(output)+1]] <- sd(depths)
  names(output) <- c("Nump", "MSI", "ED", "PS", "PD",
                     "extent", "depth", "vdepth")
  # horizontal - Total visible horizontal or terrestrial area
  # relief - Variation (Standard deviation) in elevation of the visible ground surface.
  # dsm <- terra::crop(subdsm, terra::ext(viewshed@extent, xy = TRUE))
  # dtm <- terra::crop(subdtm, terra::ext(viewshed@extent, xy = TRUE))
  dtm_z <- terra::extract(subdtm, visiblepoints)[,1]
  dsm_z <- terra::extract(subdsm, visiblepoints)[,1]
  delta_z <- dsm_z - dtm_z
  z <- cbind(dtm_z, dsm_z, delta_z)
  z <- z[which(z[,3]<=error),]
  # horizontal
  output[[length(output)+1]] <- length(z[,3]) * resolution^2
  # relief
  output[[length(output)+1]] <- sd(z[,1])
  names(output) <- c("Nump", "MSI", "ED", "PS", "PD",
                     "extent", "depth", "vdepth",
                     "horizontal", "relief")

  # horizontal visible area as SpatRaster if requested.
  if (horizontal_raster){

    mtrx <- ifelse(as.matrix((submodel<= error), wide=TRUE) &
                             (viewshed@visible == 1),
                   1, 0)
    horizontal_vis_rast <- terra::rast(
      nrows=nrow(dsm),
      ncols=ncol(dsm),
      nlyrs=1,
      extent=ext(dsm),
      vals=mtrx,
      crs=crs(dsm)
    )
    output[[length(output)+1]] <- horizontal_vis_rast
    names(output) <- c("Nump", "MSI", "ED", "PS", "PD",
                       "extent", "depth", "vdepth",
                       "horizontal", "relief", "horizontal_visible")
  }

  # skyline - Variation of (Standard deviation) of the vertical viewscape
  # (visible canopy and buildings)
  if (length(masks) == 2) {
    if (isFALSE(terra::crs(masks[[1]], proj = TRUE) == viewshed@crs)) {
      masks[[1]] <- terra::project(masks[[1]], y=terra::crs(viewshed@crs))
    }
    if (isFALSE(terra::crs(masks[[2]], proj = TRUE) == viewshed@crs)) {
      masks[[2]] <- terra::project(masks[[2]], y=terra::crs(viewshed@crs))
    }
    masks_1 <- terra::crop(masks[[1]], terra::ext(viewshed@extent, xy = TRUE))
    masks_2 <- terra::crop(masks[[2]], terra::ext(viewshed@extent, xy = TRUE))
    dsm_z <- terra::extract(dsm, visiblepoints)[,1]
    masks_1 <- terra::extract(masks_1, visiblepoints)[,1]
    masks_2 <- terra::extract(masks_2, visiblepoints)[,1]
    mask_df <- cbind(masks_1, masks_2, masks_1+masks_2, dsm_z)
    mask_ <- mask_df[which(mask_df[,3] != 0),]
    if (length(mask_df[,4]) > 1) {
      output[[length(output)+1]] <- sd(na.omit(mask_df[,4]))
    } else {
      output[[length(output)+1]] <- 0
    }
    names(output) <- c("Nump", "MSI", "ED", "PS", "PD",
                       "extent", "depth", "vdepth",
                       "horizontal", "relief", "skyline")
  }
  return(output)
}
