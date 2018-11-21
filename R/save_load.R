save_population <- function(population, file_name, compression = TRUE) {
  saveRDS(population, file = file_name, compress = compression)
}

load_population <- function(file_name) {
  readRDS(file_name)
}