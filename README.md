<p align="center">
 <img src="https://github.com/BioinfoUninaScala/epistats/blob/main/data-raw/" width="150" alt="EpiStatMOD Logo">
</p>


## EpiStatMOD
#### A new R package for the qualitative analysis of DNA methylation

### Introduction 

<p align = "justify"> 
 EpiStatMOD is a new R package aimed to extract epialleles information from any type of long- and short-read sequencing data, from targeted to genome-wide experimental data, including BS-seq, ONT-seq, Illumina 5-base seq and Biomodal 6-base seq. The epiallele-based analysis (EBA) relies on the characterization of the specific methylation patterns present on each sequenced read coming from sequencing experiments. 
 Given a genomic locus containing n Cs in the CpX context, all the possible combinations of the methylation states of these Cs observed at the single reads level are defined as epialleles. 
 Considering each sequenced read as coming from a single cell, this type of analysis can provide additional insights about the epigenetic cellular heterogeneity characterizing a sample. 
 EpiStatMOD is implemented with the ability to extract the epialleles information from genomic units that can be easily designed by the user modifying several parameters. Beyond the profiling of the epialleles composition for each analysed sample, the package is then provided with dedicated statistical functions aimed to compare epialleles compositions among different biological conditions. Finally, epialleles extraction functions are designed to easily perform strand-specific and non-CG methylation analysis. 

---------

Detailed documentation can be found at [https://bioinfouninascala.github.io/epistats](https://github.com/BioinfoUninaScala/epistats).
