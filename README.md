# Tracing 600 years of long-distance Atlantic cod trade in medieval Oslo

This repository contains the data analysis performed for the publication: *Tracing 600 years of long-distance Atlantic cod trade in medieval and post-medieval Oslo using stable isotopes and ancient DNA* by Martinez-Garcia, Pulido et al., 2024

Find the preprint of our study in [biorxiv here](https://www.biorxiv.org/content/10.1101/2024.01.25.577044v1)

The repository contanins the following files:

  A. Martinez-Garcia&Pulido.2024.Medieval_cod_trade_Oslo_ancientDNA_Isotopes.R: 
An R script that perform all statistics presented in figures and tables in the paper: *"Martinez-Garcia&Pulido.2024.Medieval_cod_trade_Oslo_ancientDNA_Isotopes.R"*
This script will read all the files in the "Input_files" folder to perform the different statistical analyses presented in the paper, as well as the plots in the main text and supplementary material.

  B. Input_files: A folder with outout files obtained from the sequencing data and the isotopic analyses of the bone material. These files are the R script presented above.

We have analysed 106 Atlantic cod bone samples collected from two archaeological sites in Oslo. Samples –cranial and postcranial elements– were morphologically (see Archive: Osteological collections, University Museum, University of Bergen) and genetically identified as Atlantic cod. Stable isotope analysis was performed on 100 Atlantic cod specimens (50 of which were also processed for genomic sequencing). The following files are the output of the isotopic as well as genomic analyses, as follow:

  1. MASTER_isotopes.genotype.txt: Results after measuring stable carbon (δ13C), nitrogen (δ15N), non-exchangeable hydrogen (δ2H, and sulphur (δ34S) on purified bone collagen. 
  2. Paleomix_summary.txt : this file contains the results produced by PALEOMIX v1.2.13, which sumarize the Clonality and Endogenous DNA fraction of the sequenced reads. Yo can also find the locality of each sample in this table.

To determine the biological origin of Atlantic cod specimens we the BAMscorer pipeline which generate probability assignments on the major chromosomal inversions in Atlantic cod (LG1, LG2, LG7 and LG12). The results and corresponding probailities are found respectively in the following files:

  3. InversionCaller_LG01.txt
  4. InversionCaller_LG02.txt
  5. InversionCaller_LG07.txt
  6. InversionCaller_LG12.txt
  7. LG01_counts.txt
  8. LG02_counts.txt
  9. LG07_counts.txt
  10. LG12_counts.txt
  11. Probability_ancient.genotypic.affinity_FoodImpact.txt

Finally, genotypes from the different inversions and their probability scores are summarized and paired with the metadata information (bone element, stratigrafic layer, Age and body size) of each individual, in the following table:

  12. All_genotypes.txt
