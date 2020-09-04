This work was performed using publicly available CAM method, and the detailed instruction of the instalment is illustrated in website:
https://bioconductor.org/packages/release/bioc/html/CAMTHC.html

To install package of CAM method that performed gene expression decomposition, start R (version "4.0") and enter:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CAMTHC")


An updated R package of CAM (renamed as deCAM) was also provided. To install this package, please enter:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("debCAM")

The detailed instruction of the instalment is illustarated in website:
http://www.bioconductor.org/packages/release/bioc/html/debCAM.html


#This is a demo program to perform gene expression data deconvolution.
#Examples of gene expression data were provided with 100 samples (Exampletrain.txt) for #training and 79 samples for testing (test.txt)
# 
# A reference gene expression matrix that represents the overall gene expression level of the subclones was estimated with an unsupervised CAM method.
#This reference matrix was further used as a guide for supervised gene expression profile deconvolution (the example data are provided as radiogenomicdata.txt)
#The genomic subclone was re-analysed using the reference matrix  

 