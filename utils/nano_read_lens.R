# Load necessary libraries
library(ShortRead)
library(ggplot2)
library(tools)

# Function to print usage information
print_usage <- function() {
  cat("Usage: Rscript nano_read_lens.R <fastq_path> <output_path>\n")
  cat("Example: Rscript nano_read_lens.R fastq/ output/\n")
}

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  print_usage()
  stop("Error: Please provide both the fastq path and output path.")
}

# Assign command line arguments to variables
fastq_path <- args[1]
output_dir <- args[2]

# List all FASTQ files directly (one file per sample)
fastq_files <- list.files(fastq_path, pattern = "\\.fastq.gz$", full.names = TRUE)

# Function to generate and save histogram for a given FASTQ file
generate_histogram <- function(fastq_file, output_file) {
  # Read the FASTQ file
  reads <- readFastq(fastq_file)
  read_lengths <- width(sread(reads))
  
  # Check if read_lengths is not empty
  if (length(read_lengths) > 0) {
    mean_length <- 10^ceiling(log10(mean(read_lengths)))
    bin_width <- mean_length / 10
    x_ticks <- seq(0, max(read_lengths), by = bin_width)
    
    # Create the histogram plot
    plot <- ggplot(data.frame(Length = read_lengths), aes(x = Length)) +
      geom_histogram(binwidth = bin_width, fill = "#ACE1AF", color = "#75A47F") +
      scale_x_continuous(breaks = x_ticks, labels = scales::comma) +
      labs(title = paste("Histogram of Read Lengths -", basename(fastq_file)), 
           x = "Read Length", y = "Frequency") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    # Save the plot as PNG
    ggsave(output_file, plot = plot, width = 10, height = 6, dpi = 300)
  } else {
    message(paste("No read lengths found in file:", fastq_file))
  }
}

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Iterate over each FASTQ file and generate the plot
for (fastq_file in fastq_files) {
  output_file <- file.path(output_dir, paste0(file_path_sans_ext(basename(fastq_file)), "_read_length_histogram.png"))
  generate_histogram(fastq_file, output_file)
}
