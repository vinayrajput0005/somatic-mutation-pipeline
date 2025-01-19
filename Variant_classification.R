# Load the necessary library
library(maftools)

# Set the working directory
setwd("/Users/vr/Desktop/Pupil_Bio/Task2/R_analysis/pupilgastk")

# Load the MAF file (Mutational Annotation Format)
maf_file <- 'PA_somatic_variants_functotated.vcf'

# Read the MAF file
tcga_luad_maf <- read.maf(maf = maf_file)

# Display the MAF data object
print(tcga_luad_maf)

# Plot the summary of the MAF file
plotmafSummary(
  maf = tcga_luad_maf,       # MAF object
  rmOutlier = TRUE,          # Remove outliers for better visualization
  addStat = 'median',        # Add median statistic
  dashboard = TRUE,          # Enable dashboard view
  titvRaw = FALSE            # Disable the raw transition/transversion plot
)
