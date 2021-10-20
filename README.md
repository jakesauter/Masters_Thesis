<!-- output:  -->
<!--   html_document:  -->
<!--     toc: true -->
<!--     toc_float: true -->

![Thesis Flowchart](figures/Thesis%20Flowchart.png)

# A Pipeline for Pan-Myeloid Flow Cytometry Data Processing and Clustering Analysis

*Presented to the Faculty of the Weill Cornell Graduate School of
Medical Sciences Cornell Universityin Partial Fulfillment of the
Requirements for the Degree of Master of Science in Computational
Biology*

by Jake Sauter May 2020

© 2021 Jake Sauter

# Table of Contents

-   [A Pipeline for Pan-Myeloid Flow Cytometry Data Processing and
    Clustering
    Analysis](#a-pipeline-for-pan-myeloid-flow-cytometry-data-processing-and-clustering-analysis)
-   [Abstract](#abstract)
-   [Biographical Sketch](#biographical-sketch)
-   [Dedication](#dedication)
-   [Acknowledgements](#acknowledgements)
-   [Contributions of Author](#contributions-of-author)
-   [List of Tables](#list-of-tables)
-   [List of Figures](#list-of-figures)
-   [List of Abbreviations](#list-of-abbreviations)
-   [Introduction](#introduction)
    -   [General Introduction](#general-introduction)
        -   [Flow Cytometry](#flow-cytometry)
        -   [Single-Cell Data Processing
            Methods](#single-cell-data-processing-methods)
    -   [Mathematical Background](#mathematical-background)
-   [Methods](#methods)
-   [Results](#results)
-   [Conclusions](#conclusions)
-   [References](#references)
-   [Code Repository](#code-repository)
-   [Notes on Process](#notes-on-process) - [Current Goal: Research
    Presentation in a
    Month](#current-goal:-research-presentation-in-a-month) - [Research
    Angle: Normal Bone Marrow
    Samples](#research-angle:-normal-bone-marrow-samples) - [Robustness
    of pipeline](#robustness-of-pipeline) - [Experimental
    Take-Aways](#experimental-take-aways) - [Notes From George’s
    Presentation](#notes-from-george's-presentation) - [Lab Meeting /
    Research Day](#lab-meeting-/-research-day) - [Possible
    Pre-Processing
    Improvements](#possible-pre-processing-improvements) - [Diffusion
    Maps](#diffusion-maps) - [Trajectory Inference /
    Pseudotime](#trajectory-inference-/-pseudotime) - [Combining flow
    initiatives](#combining-flow-initiatives) -
    [SPADE.downsampleFCS()](#spade.downsamplefcs())

# Abstract

\[Text. Indent first line of every paragraph.\]

# Biographical Sketch

Jake Sauter began his programming journey through an “Introduction to
Computer Programming in Java” class offered sophomore year at Arlington
High School in Lagrangeville New York. Jake Sauter went on to recieve
his B.S. in Applied Mathematics from State University of New York at
Oswego and is now a candidate for the Master’s of Science in
Computational Biology at Weill Cornell Graudate School of Medical
Sciences. Jake has completed his thesis work as a member of the lab of
Dr. Elli Papaemmauil in the Department of Biostatistics and Epidemeology
at Memorial Sloan Kettering Cancer Center.

# Dedication

This work is dedicated to my previous and current academic influences,
as well as to the families effected by devastating bone marrow cancers
and disorders.

# Acknowledgements

I would like to thank Elli Papaemmanuil for letting me join the lab and
work alongside the fantastic team she has put together Georgios
Asimomitis for mentoring me, always providing a wealth of ideas and
directions that could prove frutiful. I would also like to thank Ronglai
Shen for co-mentoring me and providing an outside perspective Mikhail
Roshal for providing extensive data expert pathological insight on the
Armaan Kholi for assisting with engineering tasks on the project. I am
greatful for effort you all have provided to help build my career in
Computational Biology.

# Contributions of Author

Unless specified otherwise here (and in

# List of Tables

Table 1: Title . . . . . . . . . \#\#

Table 2: Title . . . . . . . . . \#\#

# List of Figures

Figure 1: Title . . . . . . . . . \#\#

Figure 2: Title . . . . . . . . . \#\#

# List of Abbreviations

AUC Area under curve  
cfDNA Cell free DNA  
CIN Chromosomal instability  
CNV Copy number variation  
ctDNA Circulating tumor DNA  
ecDNA Extrachromosomal DNA  
MSK Memorial Sloan Kettering  
NGS Next Generation Sequencing  
ROC Receiver operating characteristic  
UMI Unique molecular index  
VAF Variant Allele Frequency  
PSA Prostate specific antigen  
NSE Neuroendocrine markers  
TF Transcription factor  
TCGA The Cancer Genome Atlas  
GES Gene expression signatures  
GSEA Gene Set Enrichment Analysis  
NES Normalized enrichment score  
ES Enrichment score  
FDR False discovery rate

# Introduction

## General Introduction

I am working on a project under MSK MIND, we have this data, I do
preprocessing, clustering and trajectory inference. The pipeline is nice
and reproducible.

### Flow Cytometry

### Single-Cell Data Processing Methods

## Mathematical Background

# Methods

# Results

# Conclusions

# References

1.  <a href="https://www.jimmunol.org/content/jimmunol/205/3/864.full.pdf" class="uri">A Comprehensive Workflow for Applying Single-Cell Clustering and Pseudotime Analysis to Flow Cytometry Data</a>

    -   Basically modelling the data-processing aspects of my thesis
        after this

2.  <a href="https://arxiv.org/abs/1305.1422" class="uri">Somoclu: An Efficient Parallel Library for Self-Organizing Maps</a>
    <https://cran.microsoft.com/snapshot/2017-09-08/web/packages/Rsomoclu/Rsomoclu.pdf>

3.  [The dynamics and regulators of cell fate decisions are revealed by
    pseudotemporal ordering of single cellsSingle-Cell Mass Cytometry of
    Differential Immune and Drug Responses Across a Human Hematopoietic
    Continuum](https://www.nature.com/articles/nbt.2859)

4.  <a href="https://www.nature.com/articles/s41596-019-0246-3.pdf?proof=t" class="uri">FLOW-MAP: a graph-based, force-directed layout algorithm for trajectory mapping in single-cell time course datasets</a>

5.  <a href="https://www.nature.com/articles/s41587-019-0071-9.pdf" class="uri">A comparison of single-cell trajectory inference methods</a>

6.  <a href="https://www.science.org/lookup/doi/10.1126/science.1198704" class="uri">Single-Cell Mass Cytometry of Differential Immune and Drug Responses Across a Human Hematopoietic Continuum</a>

7.  [Extracting a Cellular Hierarchy from High-dimensional Cytometry
    Data with SPADE cell
    assign](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3196363/)

8.  <a href="https://www.nature.com/articles/s41592-019-0529-1.pdf" class="uri">Probabilistic cell-type assignment of single-cell RNA-seq for tumor microenvironment profiling</a>

9.  [Slingshot: Cell lineage and pseudotime inference for single-cell
    transcriptomics](https://www.biorxiv.org/content/10.1101/128843v1.full)

10. [Extracting a cellular hierarchy from high-dimensional cytometry
    data with SPADE](https://pubmed.ncbi.nlm.nih.gov/21964415/)

11. Porte, J. de la, et al. “An Introduction to Diffusion Maps.” Applied
    Mathematics. Division, Department of Mathematical Sciences,
    University of Stellenbosch, South Africa,
    <a href="https://inside.mines.edu/~whereman/talks/delaPorte-Herbst-Hereman-vanderWalt-DiffusionMaps-PRASA2008.pdf" class="uri">&lt;https://inside.mines.edu/~whereman/talks/delaPorte-Herbst-Hereman-vanderWalt-DiffusionMaps-PRASA2008.pdf&gt;</a>.

12. <a href="https://pydiffmap.readthedocs.io/en/master/theory.html" class="uri">&lt;https://pydiffmap.readthedocs.io/en/master/theory.html&gt;</a>

13. <a href="https://academic.oup.com/bioinformatics/article/28/18/2400/251629" class="uri">CytoSPADE: high-performance analysis and visualization of high-dimensional cytometry data</a>

    -   Dimensionality reduction, for example, helps in reducing the
        amount of noise in the data and in visualizing the data, but a
        variety of approaches are available, with a potentially large
        impact on the final result (see Figure S1). Monocle uses
        independent component analysis, Waterfall and TSCAN use
        principal component analysis (PCA), Embedder uses Laplacian
        eigenmaps (Belkin and Niyogi, 2003), and Wishbone uses diffusion
        maps for analysis and t-distributed stochastic neighbor
        embedding (t-SNE)

14. <a href="https://www.nature.com/articles/s41587-019-0071-9.pdf" class="uri">A comparison of single-cell trajectory inference methods</a>

    -   A comparison of single-cell trajectory inference Slingshot,
        Monocle ICA vs DDRTree, PAGA, Wishbone, Wanderlust

15. <a href="https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-019-1663-x/MediaObjects/13059_2019_1663_MOESM1_ESM.pdf" class="uri">PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells - Genome Biology</a>

    1.  <a href="https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-019-1663-x/MediaObjects/13059_2019_1663_MOESM1_ESM.pdf" class="uri">PAGA Supplementary File 1 – Methods</a>

# Code Repository

Code and documentation associated with this paper can be found at the
following github repsitory:
<https://github.com/jakesauter/Masters_Thesis>

# Notes on Process

### Current Goal: Research Presentation in a Month

-   Research Day: Primary and secondary mentors are invited ⁃ Monday
    before thanksgiving ⁃ November 22nd

-   Would like to have some nice visuals that I have analyzed and some
    conclusions

### Research Angle: Normal Bone Marrow Samples

-   55 Normals ⁃ Could not find top level MRN folder for 1 patient ⁃
    More than one MRN folder found for 1 patient

### Robustness of pipeline

-   Cell ids before down-sampling?

### Experimental Take-Aways

-   SPADE diffusion maps ⁃ Behavior is not patient-specific, but method
    (SPADE) specific ⁃ Does not depend on number of neighbors

### Notes From George’s Presentation

-   MPN Diagnostic Workup ⁃ Relies on Morphology ⁃ Less than 30% blasts
    MDS vs AML

-   Flow cytometry is a useful method of measurement as it can capture
    multi-dimensional signals at a single-cell level ⁃ What is FSC, SSC
    ⁃ Conjugated antiobdies for markers of interest ⁃ M1 and M2 tube

-   Processing Pipeline ⁃ Step by step going through ⁃ Keeping table of
    content on the side and highlighting ⁃ Reference in small text at
    bottom of slide

-   Normalization ⁃ Gaussnorm not really working ⁃ Split into batches to
    only shift batch-wise and not sample-wise ⁃ Will determine shifts
    based off of “normal” samples

-   Downsampling ⁃ One option is uniform random sampling ⁃ SPADE
    density-based downsampling ⁃ Don’t over-sample cells in dense
    regions

-   Data Mining ⁃ Self Organizing Maps ⁃ Leiden Community Detection ⁃
    Diffusion Maps ⁃ UMAPs

-   Sarab Shah: ⁃ Complemented the project ⁃ Cell assign:
    <https://www.nature.com/articles/s41592-019-0529-1.pdf> ⁃ Wants us
    to talk to Toxicity project ⁃ Standardized code running to share
    with other project ⁃ Maybe run Leiden on data and subsample per
    cluster in order to subsample cells ⁃ Do you have quantification of
    cell types

-   Jacob Glass ⁃ Has been in contact with Brent Wood ⁃ Can parse blast
    from reports

### Lab Meeting / Research Day

-   Start thinking of lab meeting

-   Anybody in the project should attend

    -   Of course co-supervisor should attend

### Possible Pre-Processing Improvements

-   How did previous papers find the optimal cofactor

    -   Controls?

### Diffusion Maps

-   Can we see the eigenvalues of the diffusion components to see how
    important some are?

    -   Seems like a good piont but I think they are just comparable

    -   Maybe increases processing time because number of components is
        a paramaeter to scanpy diffusion

-   A Diffusion map is the basis transformation using the eigenvectors
    associated with the M largest eigenvalues of the Diffusion Distance
    Matrix, where the distance between coordinates of the data points in
    the new space approximates the diffusion distance, where the
    diffusion distance is the sum of the (non-normalized?) probabilites
    (density?) of the start and end points being strongly connected
    through every other point for which a path exists (p &gt; 0)

    -   Define a kernel K and calculate the kernel matrix (distances /
        non-normalized probs)

    -   Diffusion matrix = row normed K

    -   Eigenvectors of Diff Mat

    -   Map to the d-dimensional diffusion space at time t, using d
        dominant eigenvectors indicated by the magnitude of their
        eigenvalue

        -   Ex. y\_1 = a1*x\_1 + a12*x\_2 + … + a1n*x\_n y\_2 = a1*x\_1
            + a2*x\_2 + … + an*x\_n … y\_n = a1*x\_1 + a2*x\_2 + … +
            an\*x\_n

### Trajectory Inference / Pseudotime

-   Try these three methods

    -   PAGA

    -   Wishbone

    -   Wanderlust

### Combining flow initiatives

-   Sarab shah wants a flow-cyto pipeline focused meeting

-   Maybe a 30 minute presentation Wednesday 4-5pm

-   18-20 people to bring together flow initiatives

-   Table of runtime per step per sample, explain different steps

### SPADE.downsampleFCS()

-   target\_pctile — Numeric value in \[0,1\]. Densities below this
    percentile, but above exclude\_pctile will be retained. Only
    meaningful if desired\_samples is NULL.

-   desired\_samples — Desired number of samples. If set to integer
    value, the target percentile will be set internally to downsample to
    approximately the desired number of samples.

-   Seems like SPADE.downsampleFCS() only includes samples within a
    certain range of density, when it seems like in our case it would be
    best to retain all cell types, just at similar densities?

⁃ Can’t trust the Disease annotations until they are reviewed ⁃ Misha
said review genomic clusters
