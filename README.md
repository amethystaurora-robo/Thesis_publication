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

This pipeline requires input transcriptomic and metabolomic data (the data I used is confidential). Genes should be in rows and samples in columns. A separate metadata sheet was included to map column names to specific dose levels and time points. My data looked at 2 doses of ethoprophos (a high and low dose), with the transcriptome and metabolome extracted at 9 timepoints. The structure of my data is illustrated below.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/High%20dose.png"/>
</p>

A look at the overall workflow is below.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/overall_workflow.png"/>
</p>

And a look at the specific coding steps in case anyone wants to reuse the pipeline.

The first step was visualising and pre-processing the data, in files ... and ... The output from ... is then passed to the DeSeq2 analysis portion. DESeq2 is an R library used to identify differentially expressed genes (DEGs) from RNA-seq data. I used the Likelihood Ratio Test (LRT) to test the significance of the interaction between time and dose — in other words, to determine whether the effect of dose on gene expression changes over time, or vice versa.

The LRT identified 2,967 upregulated and 3,868 downregulated genes showing significant changes at the dose × time interaction level.

To interpret these results in detail, I performed Wald tests to identify DEGs at each dose level for each timepoint (a total of 18 comparisons: 9 timepoints × 2 doses). The 18 DEGs are saved to 18 files, beginning with ... I recommend saving them in a folder. The pre-processing step in ... will access the folder.

I then used WebGestalt’s Gene Set Enrichment Analysis (GSEA) to annotate these DEGs with pathway information, after pre-processing the files in ... Next, post-processing file ... is used to put all the genes in one dataset again. This file also filters the lists to the ones we care about - the ones showing significant changes at the dose x time interaction level from LRT. 

Finally, I visualised how these enriched pathways change over time and across dose levels using Tableau, combining the results into a dynamic, interpretable format. The visualisation highlights how specific biological pathways respond differently depending on the dose and the time of exposure.

<p align="left">
  <img src="https://github.com/amethystaurora-robo/Thesis_publication/blob/main/Vizzes/ES_over_time.png"/>
</p>

Note: Daphnia magna genes were mapped to pathways by aligning them to homologous human genes or KEGG Orthology (KO) terms. Due to evolutionary divergence and incomplete annotation, not all Daphnia genes could be mapped, resulting in a reduced gene set used for pathway analysis.

TODO: 
WGCNA
DynGENIE3
Fix WGCNA visual
Add versions
Add filtering at each step
Add specific file names to other specific file names
Metabolomics section
MOFA






