# Ethoprophos and Daphnia magna: A Novel Network Analysis Pipeline for the Whole Systems Modelling of a Temporal Transcriptome

## Introduction

Pesticides are used to kill pests, but they affect other organisms. As the prevalence of the organophosphate pesticide, ethoprophos, increases in the environment, risk of exposure increases for non-target organisms. Exposure to ethoprophos has been linked to physiological disturbances in earthworms, neurobehavioral impairments in rats, and in one case study, seizures in a human child.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/non-target_organisms.png">
</p>

Toxicity testing can help determine what level of exposure is "safe" in the environment. To conduct a test, a model organism (like the water flea, *Daphnia magna*, shown below) is dosed with the chemical of interest in a lab. In the past, the organism was observed and its behavioral changes, reproduction, or lifespan was recorded. Now, high-throughput 'omics technologies allow for complete extraction of gene expression data, allowing insight into cellular mechanisms which cause observed changes.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/Modern%20method.png">
</p>

However, this high-throughput data is HUGE. And many times, the data is high-dimensional, with multiple doses at multiple timepoints. A huge and complex dataset needs a complex pipeline to analyze it. In this study, I demonstrate the use of complex Graph ML alongside statistical analysis to provide a pipeline for temporal omics data. Although I won't share my full biological interpretation results here, I have proven that this pipeline is successful by backing it up with existing literature. We are working on publishing this case study!

## Methods

### Transcriptomics

This pipeline requires input transcriptomic and metabolomic data (the data I used is confidential). Genes should be in rows and samples in columns. A separate metadata sheet was included to map column names to specific dose levels and time points. My data looked at 2 doses of ethoprophos (a high and low dose), with the transcriptome and metabolome extracted at 9 timepoints. The structure of my data is illustrated below.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/High%20dose.png"/>
</p>

The transcriptomic portion of the analysis comes in 4 basic steps. 
1) Identifying significant genes
2) Annotating the genes on pathways
3) Grouping genes by their temporal co-expression patterns
4) Determining hub nodes by building a Gene Regulatory Network (GRN)

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/overall_workflow.png"/>
</p>

A more detailed look at the workflow for this data is shown below.

The specific file input and output is shown below. 

#### Step 1: Pre-processing and Visualisation
The first step was visualising and pre-processing the data. Load the data into *pre-processing_visualization.ipynb* to get an overall look at the data with violin charts, boxplots, line graphs, and PCA scatterplots. Next, raw transcriptomic data and metadata are pre-processed in *rna_preproc.r* to prepare for statistical analysis. The output of this file is *sample_sheet_filtered.csv* (metadata) and *flt_counts_filtered.csv* (read counts), used for statistical analysis. Counts are further normalized and stabilized, and the output titled *flt_counts_vst.csv* is passed to *WGCNA.r* for clustering into shared co-expression patterns.

#### Step 2: Statistical Analysis
The outputs of *rna_preproc.r* - *sample_sheet_filtered.csv* and *flt_counts_filtered.csv* are loaded into *Deseq_LRT.Rmd* for statistical testing. (Don't use the normalized and stabilized versions, since DeSeq2 has its own normalization method.)

DESeq2 is an R library used to identify differentially expressed genes (DEGs) from RNA-seq data. I used the Likelihood Ratio Test (LRT) to test the significance of the interaction between time and dose — in other words, to determine whether the effect of dose on gene expression changes over time, or vice versa. LRT fits a linear model (shown below) with and without the interaction term to determine significantly upregulated and downregulated genes at the interaction of dose and time compared to control.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/lrt.png"/>
</p>

To interpret these results in detail, I performed Wald tests to identify DEGs at each dose level for each timepoint (a total of 18 comparisons: 9 timepoints × 2 doses). The 18 DEGs are saved to 18 files, saved to a folder titled *DeSeq2 Results*. Results from LRT are saved to a file called *LRT_sig_genes.csv*. 

#### Step 3: Annotation

I then used WebGestalt’s Gene Set Enrichment Analysis (GSEA) to annotate these DEGs with pathway information. For this step, I used KEGG ontology and human orthologs for *Daphnia magna*. These were obtained following the method below.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/KOs.png"/>
</p>

*NOTE: While you obtain KOs, it will also be important to determine which are transcription factors (TFs) for the GRN building step later.*

The list of KOs obtained, and the folder of *DeSeq2 Results* will be loaded into the file *gestalt_pre-processing.ipynb*. This file is used to put the genes in the correct format for Web Gestalt. (It is also used to annotate genes for the GRN in a later step-but you can skip that for now.)

After pre-processing the files, use Web Gestalt to upload the files, selecting the correct orthologs and which information you would like. Access the website through the button below:
<p align="left">
  <a href="https://www.webgestalt.org/" target="_blank">
    <img src="https://img.shields.io/badge/Use%20WebGestalt-766090?style=for-the-badge&logo=tableau&logoColor=white"/>
  </a>
</p>

Save the results from the web tool. If uploading multiple files, the post-processing step in *gsea_processing.ipynb* will combine files into one. This file also filters the lists to the ones we care about - the ones showing significant changes at the dose x time interaction level from LRT. 

Finally, I visualised how these enriched pathways change over time and across dose levels using Tableau, combining the results into a dynamic, interpretable format. The visualisation highlights how specific biological pathways respond differently depending on the dose and the time of exposure.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/ES_over_time.png"/>
</p>

*Note: Daphnia magna genes were mapped to pathways by aligning them to homologous human genes or KEGG Orthology (KO) terms. Due to evolutionary divergence and incomplete annotation, not all Daphnia genes could be mapped, resulting in a reduced gene set used for pathway analysis.*

#### Step 4: Co-Expression Clustering
#### Step 5: (GRN)
TODO: 
WGCNA
DynGENIE3
Fix WGCNA visual
Add versions
Add filtering at each step
Add specific file names to other specific file names
Metabolomics section
MOFA

### Metabolomics






