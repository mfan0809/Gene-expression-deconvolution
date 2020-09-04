
#This is a demo program to perform gene expression data deconvolution.
#Examples of gene expression data were provided with 100 samples (Exampletrain.txt) for training and 79 samples for testing (test.txt)
# 
# A reference gene expression matrix that represents the overall gene expression level 
#of the subclones was estimated with an unsupervised CAM method.
#This reference matrix was further used as a guide for supervised gene expression profile deconvolution (the example data are provided as radiogenomicdata.txt)
#The genomic subclone was re-analysed using the reference matrix  
library(debCAM)
#read gene expression data with 100 samples
Data<-read.table("Exampletrain.txt",header= TRUE,row.names = 1)
## remove genes with zero values across 90% samples 
Data = Data[apply(Data,1,function(x) sum(x==0))<ncol(Data)*0.9,]; 
Data1<-as.matrix(Data)
set.seed(111)



rCAM <- CAM(Data1,K = 2:8,thres.low = 0.7, thres.high = 0.9)
#The MDL, a widely adopted and consistent information theoretic criterion was used to guide model selection.
#The underlying subpopulation number can be decided by minimizing the total description code length:


plot(MDL(rCAM),Data1.term=TRUE)    ##Draw minimum description length (MDL) curves at different K values
####

#############################
#############################
###The number of the subclone is selected using minimum value of MDL.
#In this example (100 samples) of gene expression data, the minimum value is achieved at 4.
minimumMDL=4;
numsubclone=minimumMDL; 
#   
#The A (proportions) and S matrix (subclone-specific expression profiles) estimated
#by CAM with a fixed subclone number

Aesttrain <- Amat(rCAM,numsubclone)    ##Estimate A matrix
MGlisttrain <- MGsforA(rCAM,K=numsubclone)    ##Obtain marker genes
simplexplot(Data1,Aesttrain,MGlisttrain)   ### Plot the simplex of the marker genes  

Sesttrain <- Smat(rCAM,numsubclone)  #Estimate S matrix


Atest<-read.table("test.txt",header= TRUE,row.names = 1)
Atest<-as.matrix(Atest)
#estimate A and S matrix from marker list 
###Molecular markers are known
Aesttest<-redoASest(Atest, MGlisttrain,S=Sesttrain,maxIter=10)$Aest

#############################################
# The trained genomic reference matrix are applied on the gene expression data 
genomicdata<-read.table("genomic2.txt",header= TRUE,row.names = 1)

genomicdata<-as.matrix(genomicdata)
#Estimating A matrix directly from known S matrix 
#The proportion matrix Aestraidogenomic was estimated using the gene expression data
#with  the Marker genes and the reference matrix from the training dataset. 
Aestradiogenomic<-redoASest(genomicdata, MGlisttrain,S=Sesttrain,maxIter=10)$Aest 

#write.table(Aestradiogenomic, file ="Aestradiogenomic.csv",sep=",") 
#write.table(MGlisttrain,file="MGlists.csv",row.names = TRUE, col.names =TRUE,sep=",")
