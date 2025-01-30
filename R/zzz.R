.onLoad <- function(libname, pkgname) {
  # Function to install and load packages
  install_and_load <- function(package) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }

  # List of required packages
  packages <- c("enrichR", "stringr", "stats", "utils")

  # Install and load required packages
  lapply(packages, install_and_load)
}
