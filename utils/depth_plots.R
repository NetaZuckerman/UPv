library(ggplot2)
library(stringr)
library(tools)

# Function to print usage information
print_usage <- function() {
  cat("Usage: Rscript depth_plots.R <depth_path> <output_path>\n")
  cat("Example: Rscript depth_plots.R depth/ output/\n")
}

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  print_usage()
  stop("Error: Please provide both the depth path and output path.")
}

# Assign command line arguments to variables
depth_path <- args[1]
output_path <- args[2]

#setwd('G:/My Drive/sheba/R_scripts/')
#depth_path <- "depth/"
#output_path <- "output/"

# Get a list of all files in the depth folder
depth_files <- list.files(path = depth_path, pattern = "\\.txt", full.names = TRUE)

# Iterate through each file
for (depth_file in depth_files) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(depth_file))
  
  # Read data from the current file
  depth <- read.csv(depth_file, sep = "\t", header = FALSE, col.names = c("ref", "Position", "Depth"))
  
  # Create the ggplot
  plot <- ggplot(depth, aes(x = Position, y = Depth)) +
    geom_area(aes(fill = Depth), color = "#F0B27A", fill = "#FAD7A0",lwd = 0.1,) +
    geom_line(y = 0, color = ifelse(depth$Depth == 0, "black", "white"), linewidth = 0.1) +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(paste(file_name,""))
  
  # Save the plot to a PNG file in the output folder
  out_file <- file.path(output_path, paste0(file_name, ".png"))
  ggsave(out_file, plot, height = 2, width = 8)
}
