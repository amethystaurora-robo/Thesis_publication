'''
This file takes the processed RNA counts and runs WGCNA to get modules.
It results in a csv with module colors ready to be visualized in Cytoscape.

'''

library(WGCNA)
library(grDevices)
library(RColorBrewer)

#load read counts and DEGs
read_counts_vst <- read.csv('rna_vst_proc.csv',header=FALSE)

head(read_counts_vst)
#process data expression for WGCNA
datExpr <- read_counts_vst[-c(1, 2), ]  # Remove the first two rows
datExpr <- datExpr[, -1]             # Remove the first column (KEGG IDs with NANs)
print(rownames(datExpr))
head(datExpr)
rownames(datExpr) <- datExpr[, 1]   # Set the second column as row names
datExpr <-datExpr[,-1] #remove first column (Dapma IDs, now rownames)
colnames(datExpr) <- datExpr[1, ]  # Set the first row as column names
head(datExpr)
datExpr <- datExpr[-1, ]            # Remove the first row
head(datExpr)
print(rownames(datExpr))
print(colnames(datExpr))
datExpr <- as.data.frame(t(datExpr))
head(datExpr,2)
sampleNames <- rownames(datExpr)
print(sampleNames)
datExpr <- datExpr [,-1]
head(datExpr,3)
length(datExpr)
print(as.vector(as.matrix(datExpr[2, ])))
#convert to matrix
datExpr <- as.matrix(datExpr)
head(datExpr)
print(rownames(datExpr))
#sample_labels <- rownames(datExpr)
length(datExpr)

#Pick soft threshold
power = pickSoftThreshold(datExpr, powerVector = seq(1, 20, by = 1), verbose = 5)

# Plotting Scale-Free topology fit index as a function of soft-thresholding power
par(mfrow = c(1, 2)) # Set up the plotting area

# Scale-Free Topology Fit
plot(power$fitIndices[, 1], power$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, R^2",
     type = "n", main = "Scale Free Topology Fit")

# Adding lines and points
text(power$fitIndices[, 1], power$fitIndices[, 2],
     labels = power$fitIndices[, 1], cex = 0.7, col = "red")

# Mean Connectivity Plot
plot(power$fitIndices[, 1], power$fitIndices[, 3],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")

# Adding lines and points
text(power$fitIndices[, 1], power$fitIndices[, 3],
     labels = power$fitIndices[, 1], cex = 0.7, col = "blue")

#Choose best soft threshold based on graphs
chosen_power <- 8

# Construct the adjacency matrix
adjacency = adjacency(datExpr, power = chosen_power)

# Create the topological overlap matrix (TOM)
#datExpr <- datExpr[,-1]
TOM = TOMsimilarity(adjacency)

#Calculate the dissimilarity TOM
dissTOM = 1 - TOM

#Create the gene tree
geneTree_degs = hclust(as.dist(dissTOM), method = "average")

#Create the dynamicMods for the subset of genes
dynamicMods_degs = cutreeDynamic(dendro = geneTree_degs, distM = dissTOM, method = "tree", minClusterSize = 50)

moduleColors <- dynamicMods_degs  #numeric module assignments
uniqueModules <- unique(moduleColors)

# Combine colors from Set1, Set2, and Set3
palette1 <- brewer.pal(9, "Set1")    # Set1 has 9 colors
palette2 <- brewer.pal(8, "Set2")    # Set2 has 8 colors
palette3 <- brewer.pal(12, "Set3")   # Set3 has 12 colors

# Combine them into one palette
combinedPalette <- c(palette1, palette2, palette3)

# Generate additional colors by interpolation if needed
combinedPalette <- colorRampPalette(combinedPalette)(35)

# Assign these colors to modules
moduleColorAssignments <- combinedPalette[match(moduleColors, uniqueModules)]

#Plot the dendrogram and module colors for the subset of genes
plotDendroAndColors(geneTree_degs, moduleColorAssignments, "Module colors", dendroLabels = FALSE, hang = 0.03)

#Print out the unique module colors and corresponding modules
print(unique(moduleColorAssignments))

# 2. Create a data frame with gene names and their corresponding module colors
geneNames = colnames(datExpr)
#sampleNames = rownames(datExpr)
print(sampleNames)
print(geneNames)
print(ncol(datExpr))
results_df = data.frame(Gene = geneNames, ModuleColor = moduleColorAssignments)

#prepare df for heatmap generation
datExpr <- datExpr[-1,]
datExpr <- apply(datExpr, 2, function(x) as.numeric(as.character(x)))
datExpr <- as.data.frame(datExpr)
MEs <- moduleEigengenes(datExpr, colors = moduleColorAssignments)$eigengenes


# Define the color names
color_names <- c("Red", "Slate Blue", "Sea Green", "Grey-Green", "Brownish Red",
                 "Dark Orange", "Bright Yellow", "Brown", "Dark Pink", "Opaque Dark Pink",
                 "Dark Blue-Green", "Turquoise", "Light Peach", "Greyish purple", "Light Purple",
                 "Beige", "Lime Green", "Opaque Yellow", "Tan", "Light Grey",
                 "Blue-Green", "Soft Green", "Washed out Yellow", "Grey-Blue", "Salmon", 
                 "Dusty Purple")

# Define the hex colors
hex_colors <- c("#E41A1C", "#556C9C", "#459D70", "#708173", "#B65C73",
                "#FF8E06", "#FFF730", "#BA7D2A", "#D56F80", "#D08AAF",
                "#8CA29B", "#6EBEA1", "#EA9369", "#AD9AAC", "#BC94C6",
                "#D0A59B", "#B5D84D", "#FFD92F", "#E9C782", "#C4B8A8",
                "#A1C2BC", "#AEDFC1", "#F7F6B7", "#C1BED7", "#EC8D8A",
                "#B29CAB")

# Create a named vector
color_map <- setNames(color_names, hex_colors)

# View the module eigengenes
head(MEs)
rownames(MEs) <- rownames(datExpr)
rownames(datExpr)

# Ensure module eigengenes and color labels are aligned
module_colors_used <- unique(moduleColorAssignments)
names(module_colors_used) <- NULL  # strip names if present

# Map hex colors to names using color_map
color_labels <- color_map[module_colors_used]

# Fallback to hex codes if no name mapping exists
color_labels[is.na(color_labels)] <- module_colors_used[is.na(color_labels)]


# Plot heatmap of Module Eigengenes
heatmap(as.matrix(MEs), 
        Rowv = NA, 
        Colv = NA, 
        scale = "row", 
        col = colorRampPalette(c("blue", "white", "red"))(50), 
        margins = c(5, 10),
        xlab = "", 
        ylab = "Module Eigengenes",
        labRow = sampleNames,
        labCol = color_labels)


write.csv(results_df,'gene_module_colors.csv',row.names=FALSE)

