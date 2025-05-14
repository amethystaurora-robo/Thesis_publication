Ethoprophos and Daphnia magna: A Novel Network Analysis Pipeline for the Whole Systems Modelling of a Temporal Transcriptome

Introduction

This study looks at the organophosphate pesticide, ethoprophos. As prevalence of ethoprophos increases in the environment, humans are increasingly at risk for exposure. You can imagine that pesticides which are used to kill organisms (in this case, the target is nematode worms) also have effects on other organisms. In this case, exposure to ethoprophos has been linked to physiological disturbances in earthworms, neurobehavioral impairments in rats, and in one case study, seizures in a human child.

"Safe" levels of chemicals are determined in part through a scientific process known as toxicity testing. A model organism (in this case the water flea, Daphnia magna) is dosed with the chemical of interest under experimental conditions. In the past, the organism was observed and its behavioral changes, reproduction, or lifespan was recorded. Since the advent of high-throughput technologies, we are now able to extract the entirety of an organism's genome or gene expression. This can give much more insight into the causes of behavioral or reproduction changes, by showing the exact cellular mechanisms which may be triggering a reaction to the chemical.

However, this high-throughput data is HUGE (in this case, over 15k genes). And many times, the data is high-dimensional, with multiple doses at multiple timepoints. A huge and complex dataset needs a complex pipeline to analyze it. This field of toxicology has been somewhat slower in adopting some of these complex pipeline formats, particularly machine learning. In this study, I demonstrate the use of complex Graph machine learning algorithms coupled with statistical analysis to provide a complete modelling pipeline for complex temporal data.

Methods

A look at the overall workflow is below.

And a look at the specific coding steps in case anyone wants to reuse the pipeline.

This pipeline requires data in a specific format. (I'm not able to share the data I used, since it belongs to an industry partner.) This is transcriptomic data from an organism, the gene names for the organism, 


