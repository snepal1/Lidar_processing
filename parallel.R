# Tree counting with lidR

library(lidR)
library(sf)
library(future)
library(future.apply)

# Filtering the noise and eliminating it
filter_noise <- function(las, sensitivity) {
  if (is(las, "LAS")) {
    p95 <- pixel_metrics(las, ~ quantile(Z, probs = 0.95), 10)
    las <- merge_spatial(las, p95, "p95")
    las <- filter_poi(las, Z > 0, Z < p95 * sensitivity)
    las$p95 <- NULL
    return(las)
  }
  if (is(las, "LAScatalog")) {
    res <- catalog_map(las, filter_noise, sensitivity = sensitivity)
    return(res)
  }
}

# Funtion to create the directory if it does not exist
create_directory_if_not_exists <- function(directory_path) {
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
    message("Created directory: ", directory_path)
  } else {
    message("Directory already exists: ", directory_path)
  }
}

# Set the number of cores to use
num_cores <- 8
# Set up parallel processing with future
plan(multisession, workers = num_cores)

cwd <- getwd()

input_data <- paste(cwd, "data", sep = "/")

# Define the LAS catalog
ctg <- readLAScatalog(input_data)

# Set the output folder for filtered LAS files
denoise <- paste(cwd, "data/denoise", sep = "/")
create_directory_if_not_exists(denoise)

# Define a function to normalize heights for a single LAS file
filter_noise_single <- function(las_file, output_folder) {
  # Normalize height using tin()
  las <- readLAS(las_file)
  las_norm <- filter_noise(las, sensitivity = 1.2)

  # Construct the output file path
  output_file <-
    file.path(output_folder, paste0("filter_", basename(las_file)))

  # Write the normalized LAS to the output folder
  writeLAS(las_norm, output_file)
}

# Use future_lapply to normalize heights in parallel
future_lapply(ctg$filename, function(las_file) {
  filter_noise_single(las_file, denoise)
})

# Define the LAS catalog
ctg <- readLAScatalog(denoise)

classify <- paste(cwd, "data/classify", sep = "/")
create_directory_if_not_exists(classify)


opt_chunk_size(ctg) <- 250 ## Chunk size 250 m, smaller the chunk faster is processing,cant go below 250
opt_chunk_buffer(ctg) <- 20 ## Buffer size optional
opt_output_files(ctg) <- file.path(classify,"Ground_{i}") #output folder where you want to keep your classified graound data, .las format
classified_ctg <- classify_ground(ctg, csf()) # classify the ground





# Define a function to normalize heights for a single LAS file
#classify_ground_height_single <- function(las_file, output_folder) {
  # Classify ground using csf
  #las <- readLAS(las_file)
  #las_norm <- classify_ground(las, csf())

  # Construct the output file path
  #output_file <-
    #file.path(output_folder, paste0("classified_", basename(las_file)))

  # Write the normalized LAS to the output folder
  #writeLAS(las_norm, output_file)
#}

# Use future_lapply to normalize heights in parallel
#future_lapply(ctg$filename, function(las_file) {
  #classify_ground_height_single(las_file, classify)
#})

# Define the LAS catalog
ctg <- readLAScatalog(classify)

# Set the output folder for normalized LAS files
normalize <- paste(cwd, "data/normalize", sep = "/")
create_directory_if_not_exists(normalize)

# Define a function to normalize heights for a single LAS file
normalize_height_single <- function(las_file, output_folder) {
  # Normalize height using tin()
  las <- readLAS(las_file)
  las_norm <- normalize_height(las, tin())

  # Construct the output file path
  output_file <-
    file.path(output_folder, paste0("norm_", basename(las_file)))

  # Write the normalized LAS to the output folder
  writeLAS(las_norm, output_file)
}

# Set the future seed option to TRUE
options(future.seed = TRUE)
# Use future_lapply to normalize heights in parallel
future_lapply(ctg$filename, function(las_file) {
  normalize_height_single(las_file, normalize)
})

ctg <- readLAScatalog(normalize)
segment <- paste(cwd, "data/segment", sep = "/")
create_directory_if_not_exists(segment)

# Define a function to normalize heights for a single LAS file
process_lidar_file <- function(las_file, output_folder) {
  # Read the LAS file
  las <- readLAS(las_file)

  # Perform processing without creating large objects in the global environment
  chm <- rasterize_canopy(las, res = 1, algorithm = p2r(subcircle = 0.2))
  ttops <- locate_trees(las, lmf(ws = 15))
  algo <- dalponte2016(chm, ttops)
  las <- segment_trees(las, algo)
  crowns <- crown_metrics(las, func = .stdmetrics, geom = "convex")

  # Extract the index from the LAS file name
  las_index <- sub(".las", "", basename(las_file))

  # Construct the output file path with the index
  output_file <- file.path(output_folder, paste0("Tree_", las_index, ".kml"))

  st_write(crowns, output_file, row.names = FALSE)
}

# Use future_lapply to process LAS files and save results in parallel
future_lapply(ctg$filename, function(las_file) {
  process_lidar_file(las_file, segment)
})
