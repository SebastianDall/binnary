# List of required packages
required_packages <- c("optparse", "tidyverse", "seqinr")  # Add all necessary packages here

# Function to install packages
install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
    }
  }
}

# Install the required packages
install_packages(required_packages)