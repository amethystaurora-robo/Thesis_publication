Ethoprophos and Daphnia magna: A Novel Network Analysis Pipeline for the Whole Systems Modelling of a Temporal Transcriptome

Introduction

This study looks at the organophosphate pesticide, ethoprophos. As prevalence of ethoprophos increases in the environment, humans are increasingly at risk for exposure. You can imagine that pesticides which are used to kill organisms (in this case, the target is nematode worms) also have effects on other organisms. In this case, exposure to ethoprophos has been linked to physiological disturbances in earthworms, neurobehavioral impairments in rats, and in one case study, seizures in a human child.

"Safe" levels of chemicals are determined in part through a scientific process known as toxicity testing. A model organism (in this case the water flea, Daphnia magna) is dosed with the chemical of interest under experimental conditions. In the past, the organism was observed and its behavioral changes, reproduction, or lifespan was recorded. Since the advent of high-throughput technologies, we are now able to extract the entirety of an organism's genome or gene expression. This can give much more insight into the causes of behavioral or reproduction changes, by showing the exact cellular mechanisms which may be triggering a reaction to the chemical.

However, this high-throughput data is HUGE (in this case, over 15k genes). And many times, the data is high-dimensional, with multiple doses at multiple timepoints. A huge and complex dataset needs a complex pipeline to analyze it. This field of toxicology has been somewhat slower in adopting some of these complex pipeline formats, particularly machine learning. In this study, I demonstrate the use of complex Graph machine learning algorithms coupled with statistical analysis to provide a complete modelling pipeline for complex temporal data.

Methods

A look at the overall workflow is below.

And a look at the specific coding steps in case anyone wants to reuse the pipeline.

This pipeline requires data input data-I'm not allowed to share the data that I used. It makes use of transcriptomic and metabolomic data, although the majority of the analysis is conducted on the transcriptome. The transcriptomic data included 29134 genes, with genes in rows and samples in columns. The samples in this case were a combination of dose level of ethoprophos (a high dose ethoprophos and a low dose ethoprophos), dosed on the organism at 9 timepoints. A separate metadata sheet was included to map column names to specific dose levels and time points. The structure of my data is illustrated below.

The first step was visualising and pre-processing the data, in files ... and ... The output from ... is then passed to the DeSeq2 analysis portion. DESeq2 is an R library used to identify differentially expressed genes (DEGs) from RNA-seq data. I used the Likelihood Ratio Test (LRT) to test the significance of the interaction between time and dose — in other words, to determine whether the effect of dose on gene expression changes over time, or vice versa.

The LRT identified 2,967 upregulated and 3,868 downregulated genes showing significant changes at the dose × time interaction level.

To interpret these results in detail, I performed Wald tests to identify DEGs at each dose level for each timepoint (a total of 18 comparisons: 9 timepoints × 2 doses).

I then used WebGestalt’s Gene Set Enrichment Analysis (GSEA) to annotate these DEGs with pathway information.

Finally, I visualised how these enriched pathways change over time and across dose levels using Tableau, combining the results into a dynamic, interpretable format. The visualisation highlights how specific biological pathways respond differently depending on the dose and the time of exposure.

Note: Daphnia magna genes were mapped to pathways by aligning them to homologous human genes or KEGG Orthology (KO) terms. Due to evolutionary divergence and incomplete annotation, not all Daphnia genes could be mapped, resulting in a reduced gene set used for pathway analysis.

TODO: 
WGCNA
DynGENIE3
Fix WGCNA visual
Add DeSeq2 visuals
Add versions
Add filtering at each step
Add specific file names to other specific file names
Metabolomics section
MOFA






