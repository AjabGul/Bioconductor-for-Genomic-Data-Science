# Install AnnotationHub and dependencies
# BiocManager::install("AnnotationHub")


# Load the AnnotationHub package to access genomic annotation datasets
library(AnnotationHub)

 

# Create an AnnotationHub object 'ah' and connect to the hub snapshot
ah <- AnnotationHub()  

# Print a summary of all available datasets in AnnotationHub
ah  

# Access the first dataset in the hub (useful to inspect a single record)
ah[1]  

# List all unique data providers available in the hub
unique(ah$dataprovider)  

# List all unique species available in the hub
unique(ah$species)  

# Subset the hub to only include datasets for human (Homo sapiens)
ah <- subset(ah, species == "Homo sapiens")  

# Print the subset to check only human datasets remain
ah  

# Query the hub to find datasets related to H3K4me3 marks in the Gm12878 cell line
query(ah, c("H3K4me3", "Gm12878"))  

# Open the AnnotationHub object in RStudio's spreadsheet-like viewer for browsing
View(ah)  


