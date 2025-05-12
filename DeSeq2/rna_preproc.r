# Loading the Libraries and Setting the Working Path
# ---------------------------------------------------
rm(list = ls())

library(ggplot2)
library(plotly)
library(reshape2)
library(sva)
library(DESeq2)



# Loading the Sample Sheet and Read Counts
# ---------------------------------------------------
# Load sample sheet
sample_sheet_file <- 'daphnia_metadata_sanitized.tsv'
sample_sheet <- read.table(sample_sheet_file, sep = '\t', header = T, row.names = 20, stringsAsFactors = T, check.names = F, comment.char = '')

sample_sheet <- sample_sheet[!(sample_sheet$compound %in% c('EXTRACTION BLANK','EMPTY_WELL')),]
sample_sheet <- droplevels(sample_sheet)

sample_sheet$timepoint <- relevel(sample_sheet$timepoint, ref = 'TP1')
sample_sheet$dose_code <- relevel(sample_sheet$dose_code, ref = 'Z') 

head(sample_sheet)
dim(sample_sheet)


# Load read counts
# raw_counts_file <- 'htseq_raw_counts.csv'
raw_counts_file <- 'kallisto_raw_counts.csv'
raw_counts <- read.table(raw_counts_file, sep = ',', header = T, row.names = 1, check.names = F, comment.char = '')

head(raw_counts)
dim(raw_counts)

# Convert to a matrix (and round to integer)
raw_counts <- as.matrix(raw_counts)
raw_counts <- round(raw_counts)


# Match the order of sample sheet and read counts
shared_ids <- intersect(colnames(raw_counts), row.names(sample_sheet))
length(shared_ids)

raw_counts <- raw_counts[,shared_ids]
sample_sheet <- sample_sheet[shared_ids,]

stopifnot(length(shared_ids) > 0)
stopifnot(row.names(sample_sheet) == colnames(raw_counts))
stopifnot(!any(is.na(row.names(sample_sheet))))
stopifnot(!any(is.na(colnames(raw_counts))))

dim(raw_counts)



# Checking the Data Quality
# ---------------------------------------------------
# plot the counts distribution of each sample
# colnames(raw_counts)
# p <- ggplot(as.data.frame(log2(raw_counts+1)), aes(x = XL100)) + 
#   geom_histogram(fill = '#619CFF', binwidth = 0.5)
# ggplotly(p)
# 
# df = melt(log2(raw_counts + 1), varnames = c('genes', 'samples'))
# df = data.frame(df, condition = sample_sheet[df$samples,]$compound)
# 
# p <- ggplot(df, aes(x = value, colour = samples, fill = samples)) +
#   geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ condition) +
#   theme(legend.position = 'top') + xlab('log2(count + 1)')
# 
# p <- ggplot(df, aes(x = samples, y = value, fill = condition)) + 
#   geom_boxplot() + 
#   ylab('log2(count + 1)')
# p
# ggplotly(p)


# Counts and null counts
df <- data.frame(
  # sample = sample_sheet$label,
  sample = row.names(sample_sheet),
  num.counts = colSums(raw_counts), 
  null.counts = colSums(raw_counts==0)/dim(raw_counts)[1]*100, 
  treatment = sample_sheet$compound
)

p <- ggplot(df, aes(x = sample, y = num.counts, fill = treatment)) +
  geom_bar(stat = 'identity') + ylab('total number of read counts') + #scale_fill_manual(values = c('#619CFF', '#F564E3')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplotly(p)

p <- ggplot(df, aes(x = sample, y = null.counts, fill = treatment)) +
  geom_bar(stat = 'identity') + ylab('percentage of null counts (%)') + #scale_fill_manual(values = c('#619CFF', '#F564E3')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplotly(p)


# plot the counts distribution of each sample
# df = melt(log2(raw_counts + 1), varnames = c('genes', 'samples'))
# df = data.frame(df, condition = sample_sheet[df$samples,]$compound)
# p <- ggplot(df, aes(x = value, colour = samples, fill = samples)) +
#   geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ condition) +
#   theme(legend.position = 'top') + xlab('log2(count + 1)')
# ggplotly(p)

# p <- ggplot(as.data.frame(log2(raw_counts+1)), aes(x = DAA001MB1)) +
#   geom_histogram(fill = '#619CFF', binwidth = 0.5)
# ggplotly(p)


# Filtering low quality samples
rm.samples <- row.names(df[df$null.counts > 75,])
# rm.samples <- row.names(df[df$num.counts < 1E+6 | df$null.counts > 80,])
# rm.samples <- c('DAA001MB4', rm.samples) #DAA006LC4
rm.samples
length(rm.samples)

raw_counts <- raw_counts[,!(row.names(sample_sheet) %in% rm.samples)]
sample_sheet <- sample_sheet[!(row.names(sample_sheet) %in% rm.samples),]
stopifnot(row.names(sample_sheet) == colnames(raw_counts))

dim(raw_counts)


# Remove low / zero variance genes
read_vars <- rowVars(raw_counts, useNames = F)
raw_counts <- raw_counts[which(read_vars > 1e-10),]
dim(raw_counts)


# Remove genes with too small total counts
read_cnts <- rowSums(raw_counts)
raw_counts <- raw_counts[read_cnts > 10,]
dim(raw_counts)


# Remove genes with too many null counts
cpm_valids <- rowSums(raw_counts/colSums(raw_counts)*1e+6 > 0.5)
raw_counts <- raw_counts[cpm_valids > 55,] # at least one pass per treatment group
dim(raw_counts)


# export
flt_sample_sheet_file <- gsub('\\..*', '_filtered.csv', sample_sheet_file)
write.csv(sample_sheet, file = flt_sample_sheet_file)

flt_counts_file <- gsub('\\..*', '_filtered.csv', raw_counts_file)
write.csv(raw_counts, file = flt_counts_file)



# Batch correction if needed
# ---------------------------------------------------
# batch <- as.integer(sample_sheet$exposure_batch)
# grp1 <- sample_sheet$compound
# grp2 <- sample_sheet$timepoint
# grp3 <- sample_sheet$dose_code
# grps <- cbind(grp1,grp2,grp3)
# raw_counts <- ComBat_seq(raw_counts, batch = batch, group = NULL, covar_mod = grps)



# Applying DESeq2 normalisation and transforrmation without fitting model
# ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
    raw_counts, colData = sample_sheet,
    design = ~ timepoint * dose_code
)


# Obtain normalised read counts
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = T)

norm_counts <- norm_counts[which(rowVars(norm_counts, useNames = F) > 1e-10),]
dim(norm_counts)

norm_counts_file <- gsub('\\..*', '_norm.csv', raw_counts_file)
write.csv(norm_counts, file = norm_counts_file)


# Obtain vst transformed read counts
vds <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_counts <- as.matrix(assay(vds))

vst_counts <- vst_counts[which(rowVars(vst_counts, useNames = F) > 1e-10),]
dim(vst_counts)

vst_counts_file <- gsub('\\..*', '_vst.csv', raw_counts_file)
write.csv(vst_counts,  file = vst_counts_file)

