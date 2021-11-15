![Thesis Flowchart](figures/Thesis_Flowchart.png)

# A Pipeline for Pan-Myeloid Flow Cytometry Data Processing and Clustering Analysis

*Presented to the Faculty of the Weill Cornell Graduate School of
Medical Sciences Cornell Universityin Partial Fulfillment of the
Requirements for the Degree of Master of Science in Computational
Biology*

by Jake Sauter May 2020

© 2021 Jake Sauter

# Table of Contents

    ## ℹ Sourcing https://gist.githubusercontent.com/gadenbuie/c83e078bf8c81b035e32c3fc0cf04ee8/raw/57e7c1a8caed373d9cca7e248cb972d01d26678a/render_toc.R

    ## ℹ SHA-1 hash of file is d520b023f6a9715c381801b34fb1942629f2df12

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
    -   [Myeloid Disorders](#myeloid-disorders)
    -   [Flow Cytometry](#flow-cytometry)
    -   [Clustering](#clustering)
-   [Methods](#methods)
    -   [Flow Cytometry Data
        Processing](#flow-cytometry-data-processing)
    -   [Clustering](#clustering)
    -   [Dimensionality Reduction](#dimensionality-reduction)
    -   [Trajectory Inference](#trajectory-inference)
-   [Results](#results)
-   [Conclusions](#conclusions)
-   [References](#references)
-   [Code Repository](#code-repository)
-   [Notes](#notes)
    -   [Research Presentation in a
        Month](#research-presentation-in-a-month)
    -   [Normal Bone Marrow Samples](#normal-bone-marrow-samples)
    -   [Robustness of pipeline](#robustness-of-pipeline)
    -   [Experimental Take-Aways](#experimental-take-aways)
    -   [Notes From MIND Presentation](#notes-from-mind-presentation)
    -   [Lab Meeting / Research Day](#lab-meeting--research-day)
    -   [Possible Pre-Processing
        Improvements](#possible-pre-processing-improvements)
    -   [Diffusion Maps](#diffusion-maps)
    -   [Trajectory Inference /
        Pseudotime](#trajectory-inference-/-pseudotime)
    -   [Combining flow initiatives](#combining-flow-initiatives)
    -   [SPADE.downsampleFCS()](#spade.downsamplefcs())
    -   [Observations in Cohort](#observations-in-cohort)
    -   [Normal Samples Diffusion Map Colored by
        Leiden](#normal-samples-diffusion-map-colored-by-leiden)
    -   [Normal Samples PAGA graph](#normal-samples-paga-graph)
    -   [Research Day](#research-day)
        -   [Next Steps](#next-steps)
        -   [Prioritization](#prioritization)
        -   [](#)
        -   [Wishbone (Wanderlust++)](#wishbone-(wanderlust++))
        -   [Wanderlust](#wanderlust)
        -   [Wishbone](#wishbone)
    -   [Monday November 11th](#monday-november-11th)
        -   [](#)
        -   [Optimal Trajectory Inference Methods for our
            case](#optimal-trajectory-inference-methods-for-our-case)
    -   [](#)
    -   [Research Day Presentation](#research-day-presentation)
        -   [FlowJo License](#flowjo-license)

# Abstract

\[Text. Indent first line of every paragraph.\]

# Biographical Sketch

Jake Sauter recieved his B.S. in Applied Mathematics from State
University of New York at Oswego and is now a candidate for the Master’s
of Science in Computational Biology at Weill Cornell Graudate School of
Medical Sciences. Jake has completed his thesis work as a member of the
lab of Dr. Elli Papaemmauil in the Department of Biostatistics and
Epidemeology at Memorial Sloan Kettering Cancer Center.

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

I am working on a project under MSK MIND, we have this data, I do
preprocessing, clustering and trajectory inference. The pipeline is nice
and reproducible.

## Myeloid Disorders

## Flow Cytometry

## Clustering

# Methods

## Flow Cytometry Data Processing

Flow cytometry is a method to study …

## Clustering

Clustering is an unsupervised machine learning technique …

## Dimensionality Reduction

Dimensionality reduction is the method of projecting / visualizing a
higher dimensional space in a lower dimensional space such that the
maximal amount of information (possibly being relations between data
points) is retained

Dimensionality reduction could prove fruitful in the field of flow
cytometry and pathology as pathologists determine cell populations and
disease state through a series of one or two-dimensonal plots in which
multi-variable patterns could be missed

## Trajectory Inference

Trjaectory inference can be used independently or alongside clustering
analysis to

# Results

# Conclusions

# References

1.  <a href="https://www.jimmunol.org/content/jimmunol/205/3/864.full.pdf" class="uri">A Comprehensive Workflow for Applying Single-Cell Clustering and Pseudotime Analysis to Flow Cytometry Data</a>

    -   Helpful blueprint for flow cytometry data-processing

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

8.  [Diffusion maps for high-dimensional single-cell analysis of
    differentiation
    data](https://academic.oup.com/bioinformatics/article/31/18/2989/241305)

9.  <a href="https://www.nature.com/articles/s41592-019-0529-1.pdf" class="uri">Probabilistic cell-type assignment of single-cell RNA-seq for tumor microenvironment profiling</a>

10. [Slingshot: Cell lineage and pseudotime inference for single-cell
    transcriptomics](https://www.biorxiv.org/content/10.1101/128843v1.full)

11. [Extracting a cellular hierarchy from high-dimensional cytometry
    data with SPADE](https://pubmed.ncbi.nlm.nih.gov/21964415/)

12. <a href="https://inside.mines.edu/~whereman/talks/delaPorte-Herbst-Hereman-vanderWalt-DiffusionMaps-PRASA2008.pdf" class="uri">An Introduction to Diffusion Maps</a>.

13. <a href="https://pydiffmap.readthedocs.io/en/master/theory.html" class="uri">pydiffmap Diffusion Maps Theory</a>

14. <a href="https://academic.oup.com/bioinformatics/article/28/18/2400/251629" class="uri">CytoSPADE: high-performance analysis and visualization of high-dimensional cytometry data</a>

15. <a href="https://www.nature.com/articles/s41587-019-0071-9.pdf" class="uri">A comparison of single-cell trajectory inference methods</a>

16. <a href="https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-019-1663-x/MediaObjects/13059_2019_1663_MOESM1_ESM.pdf" class="uri">PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells - Genome Biology</a>

    1.  <a href="https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-019-1663-x/MediaObjects/13059_2019_1663_MOESM1_ESM.pdf" class="uri">PAGA Supplementary File 1 – Methods</a>

17. <a href="https://onlinelibrary.wiley.com/doi/10.1002/eji.201646632" class="uri">Guidelines for the use of flow cytometry and cell sorting in immunological studies*</a>

18. <a href="https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.23664" class="uri">Implementation and Validation of an Automated Flow Cytometry Analysis Pipeline for Human Immune Profiling</a>

19. <a href="https://onlinelibrary.wiley.com/doi/10.1002/cyto.b.20554" class="uri">Elucidation of seventeen human peripheral blood B-cell subsets and quantification of the tetanus response using a density-based method for the automated identification of cell populations in multidimensional flow cytometry data</a>

20. [Exploration of Cell Development Pathways through High-Dimensional
    Single Cell Analysis in Trajectory
    Space](https://www.cell.com/iscience/pdf/S2589-0042(20)30025-0.pdf)

21. [Algebraic approach to single-pushout graph
    transformatio](https://www.researchgate.net/publication/223806275_Algebraic_approach_to_single-pushout_graph_transformation)

