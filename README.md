# Analysis of Madingley Food webs

This repository contains the code and data that can be used to reproduce the results of the following publication:

Cesar O. Flores, Susanne Kortsch, Derek Tittensor, Mike Harfoot, and Drew Purves.
**Food Webs: Insights from a General Ecosystem Models (2019).** 
*bioRxiv*.

For more information about the Madingley Model see: https://madingley.github.io/


## Overview

The **code** is composed of:
* **main_code:** A a collection of R functions for the analsyis of Madingley food web data. Among the functionality included
are:
  * Aggregation of Madingley output into a smaller number of nodes using K-means and Hierarchical clustering.
  * Weightening of edges according to different metrics, with the biomass transfer as the one using in the publication. 
  * Sampling of edges using three different methods, with the Poisson probabilities method as the one used in the
    publication.
* **scripts_code:** A set of R notebooks that can be used to reproduce all plots and results described in the publication.

The **data** is composed of:
* **madingley_data**: A ziped file that include node and edge data across different time snapshots produced by the Madingley Model.
* **empirial_data**: A set of empirical food web datasets.   
  

## Instructions

Before trying to run the R notebooks located under **code/scripts_code**, the following steps are required:
 * Clone the this repository in a local directory.
 * Unzip the contents of **data/madingley_data/need_to_unzip_here.zip** on **data/madingley_data/**
 * Install the following R libraries:
   * igraph
   * data.table
   * ggplot2
   * gridExtra
   * NetIndices
   * parallel (optional and only required if experiments need to run in parallel)
   * latex2exp 
 
After completing the previous steps, the user is ready to run the R notebooks (Rmd files) to reproduce our results.
Furthermore, these notebooks can be seen as examples of how to use the main code to analyze Madingley food web data. 
    