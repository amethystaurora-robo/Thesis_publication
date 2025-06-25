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

#### Step 1: Pre-processing and Visualisation
An overall look at the file order for pre-processing, statistical analysis and annotation is shown below:
<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/pre-proc-viz-methods.png"/>
</p>

The first step was visualising and pre-processing the data. Load the data into *pre-processing_visualization.ipynb* to get an overall look at the data with violin charts, boxplots, line graphs, and PCA scatterplots. Next, raw transcriptomic data and metadata are pre-processed in *rna_preproc.r* to prepare for statistical analysis. The output of this file is *sample_sheet_filtered.csv* (metadata) and *flt_counts_filtered.csv* (read counts), used for statistical analysis. Counts are further normalized and stabilized, and the output titled *rna_vst_proc.csv* is passed to *WGCNA.r* for clustering into shared co-expression patterns.

#### Step 2: Statistical Analysis
The outputs of *rna_preproc.r* - *sample_sheet_filtered.csv* and *flt_counts_filtered.csv* are loaded into *Deseq_LRT.Rmd* for statistical testing. (Don't use the normalized and stabilized versions, since DeSeq2 has its own normalization method.)

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_stat_analysis.png"/>
</p>

DESeq2 is an R library used to identify differentially expressed genes (DEGs) from RNA-seq data. I used DESeq2 version 1.38.3. I used the Likelihood Ratio Test (LRT) to test the significance of the interaction between time and dose — in other words, to determine whether the effect of dose on gene expression changes over time, or vice versa. LRT fits a linear model (shown below) with and without the interaction term to determine significantly upregulated and downregulated genes at the interaction of dose and time compared to control.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/lrt.png"/>
</p>

To interpret these results in detail, I performed Wald tests to identify DEGs at each dose level for each timepoint (a total of 18 comparisons: 9 timepoints × 2 doses). The 18 DEGs are saved to 18 files, saved to a folder titled *DeSeq2 Results*. Results from LRT are saved to a file called *LRT_sig_genes.csv*. 

#### Step 3: Annotation

I then used WebGestalt’s Gene Set Enrichment Analysis (GSEA) to annotate these DEGs with pathway information. I used Web Gestalt version 2024. For this step, I used KEGG ontology and human orthologs for *Daphnia magna*. These were obtained following the method below.

<p align="left" style="font-size: 16px;">
  1. List of gene names
  &nbsp;→&nbsp;
  <a href="http://wfleabase.org/" target="_blank">
    <img src="https://img.shields.io/badge/wFleaBase-018575?style=for-the-badge&logo=biomart&logoColor=white"/>
  </a>
  &nbsp;→&nbsp;
  3. Amino acid sequences
  &nbsp;→&nbsp;
  <a href="https://www.kegg.jp/blastkoala/" target="_blank">
    <img src="https://img.shields.io/badge/BlastKOALA%20(KEGG)-018575?style=for-the-badge&logo=tableau&logoColor=white"/>
  </a>
  &nbsp;→&nbsp;
  5. KOs
</p>

*NOTE: While you obtain KOs, it will also be important to determine which are transcription factors (TFs) for the GRN building step later.*

The list of KOs obtained, and the folder of *DeSeq2 Results* will be loaded into the file *gestalt_pre-processing.ipynb*. This file is used to put the genes in the correct format for Web Gestalt. (It is also used to annotate genes for the GRN in a later step-but you can skip that for now.)

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_annotation.png"/>
</p>

After pre-processing the files, use Web Gestalt to upload the files, selecting the correct orthologs and which information you would like. Access the website through the button below:
<p align="left">
  <a href="https://www.webgestalt.org/" target="_blank">
    <img src="https://img.shields.io/badge/Use%20WebGestalt-018575?style=for-the-badge&logo=tableau&logoColor=white"/>
  </a>
</p>

Save the results from the web tool. If uploading multiple files, the post-processing step in *gsea_processing.ipynb* will combine files into one. This file also filters the lists to the ones we care about - the ones showing significant changes at the dose x time interaction level from LRT. 

Finally, I visualised how these enriched pathways change over time and across dose levels using Tableau, combining the results into a dynamic, interpretable format. The visualisation highlights how specific biological pathways respond differently depending on the dose and the time of exposure.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/ES_over_time.png"/>
</p>

*Note: Daphnia magna genes were mapped to pathways by aligning them to homologous human genes or KEGG Orthology (KO) terms. Due to evolutionary divergence and incomplete annotation, not all Daphnia genes could be mapped, resulting in a reduced gene set used for pathway analysis.*

#### Step 4: Co-Expression Clustering
An important aspect of the temporal analysis was grouping genes by their shared co-expression patterns. In other words, these groups of genes are upregulated or downregulated together. This method uses a hierarchical clustering algorithm to group genes into 'modules.'

Normalized and vst transformed genes, in the file *rna_vst_proc.csv* is passed to *WGCNA.r* to assign genes to different modules. This produces an output file, *genes_module_colors.csv* which is passed to *network_processing.ipynb* together with outputs from the GRN for a final visualisation. I used WGCNA version 1.72.517 to perform this analysis.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_WGCNA.png"/>
</p>

In *WGCNA.R*, a heatmap is also generated to visualise shared co-expression patterns (eigengenes) by modules.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/module_eigengenes.png"/>
</p>

#### Step 5a: GRN Pre-processing
The core objective of this pipeline is to construct a Gene Regulatory Network (GRN), which maps regulatory relationships between genes. In this network, transcription factors (TFs) are key regulators that activate or repress target genes. Identifying TFs as hub nodes—those with a high number of connections—can highlight potential pathway initiators. Such TFs may serve as biomarkers of ethoprophos exposure, providing insights into its mechanism of action at the molecular level.

Because the data is temporal, a dynamic gene network inference method was used. DynGENIE3 works with temporal data to build a GRN. 

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/GRN_pre-processing.png"/>
</p>

TFs obtained in Step 3 and filtered and normalized transcriptomic data in *rna_vst_proc.csv* are loaded into *GRN_pre-processing.ipynb*. This will put the files into the correct format based on dose levels, and output necessary files including the gene order, in *gene_order.csv*, the list of known TFs in *regulatory_genes.csv*, and the gene expression levels, split into their conditions, in this case titled *control_df.txt*, *low_df.txt* and *high_df.txt*. 

#### Step 5b: Running the GRN

The wrapper doc available from https://github.com/vahuynh/dynGENIE3 will need to be compiled for use in the R file, *dynGENIE3.R*. Weights are saved in files starting with *link_list* and decay rates are stored in files starting with *alphas.* These files are loaded into *DynGENIE3_analysis.R* for filtering, and finally combined in *process_results.R*. 

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_GRN1.png?raw=true"/><img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_GRN2.png?raw=true"/>
</p>

The algorithm DynGENIE3 uses ensembles of regression trees are used to determine strengths of regulatory relationship between each TF and target gene. This model results in a list of weights - the relationship strength, and messenger RNA (mRNA) decay rates. Decay rates of mRNA are inferred by applying an exponential decay equation to the time series data. Faster decay rates lead to repression of mRNA transcription, affecting the level of proteins produced and ultimately the pathways which are enriched from these genes. Focusing on higher decay rates highlights genes which may be responding to cellular stress in this way.

DynGENIE3 has two parameters which may be tuned-k and ntrees. I have provided a file, *parameter_tuning.R* to optimize the algorithm for any data. An HPC would be beneficial for running DynGENIE3, but to perform parameter tuning, it will be completely necessary. Additionally, the file *paramter_tuning_visualisation.ipynb* can be used to create the pretty 3D mesh plot below.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/3dmesh_plot.png"/>
</p>

#### Step 5c: GRN annotation and post-processing

I used Web Gestalt to annotate the results from the GRN with pathway information, in the same process as Step 3. In *network_processing.ipynb*, the annotated genes are merged with modules in *gene_module_colors.csv* from *WGCNA.R*. 

<p align="center">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_GRN_Annotation1.png?raw=true" width="49%" /><img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/workflow_GRN_annotation2.png?raw=true" width="49%" />
</p>

The merged file *tfs_wgcna.csv* can be loaded into Cytoscape software. Below are two networks I created, annotating with pathways and transcription factor identifiers. The colors represent the co-expression groupings. Cytoscape provides its own network inference methods, by which you can find 'hub nodes'. This will allow for biological interpretation of the results. 

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/es_network_img.png"/>
</p>

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/alphas_network_img.png"/>
</p>

### Metabolomics

For this pipeline, the main focus has been transcriptomic data. However, metabolomic data is still analysed using statistical tests. Both omics levels are connected through the use of multi-omics factor analysis (MOFA). 

#### Step 1: Statistical Analysis

Metabolomic data with metabolites in rows and samples in columns is cleaned, dropping the header and columns with ortholog 
information to prepare for statistical analysis.

This script performs one-way ANOVA and post-hoc Tukey's HSD tests on positive and negative datasets
(pos_melted_renamed and neg_melted_renamed) to assess whether there are statistically significant 
differences in mass measurements across three categorical variables: Condition, Treatment, and Timepoint.

Steps:
1. Group data by each categorical variable for both positive and negative datasets.
2. Perform one-way ANOVA to test for overall group differences.
3. If significant, perform Tukey's HSD to identify specific group differences.

Outputs:
- F-statistic and p-values from ANOVA
- Tukey's test results for all pairwise comparisons

#### Step 2: Combined Analysis (multi-omics)

Multi-omics factor analysis is a statistical framework used to integrate omics datasets - in this case, transcriptomic and metabolomic datasets - by identifying latent factors that explain shared variation across datasets. This works to:

- Reduce dimensionality across multiple omics layers.
- Identify shared and dataset-specific sources of variation.
- Help to uncover biological signals, sample subgroups, or batch effects.

<p align="left">
<img src ="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/mofa_variance_explained_bar.png">
<img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/mofa_variance_explained_scatter.png">
<img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/mofa_factors.png">
</p>



