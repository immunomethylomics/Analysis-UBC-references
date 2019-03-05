#####################################################
# Author: Lucas A Salas
# Date: 06/26/2018
# The IDOL objects have been saved for reproducibilit in github and zenodo
#####################################################

#Reading the British Columbia data
setwd("D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive")
library(minfi)
library(FlowSorted.Blood.EPIC)
library(ENmix)

sheet <- read.metharray.sheet(getwd())
megdata <- read.metharray.exp(targets = sheet,extended = TRUE)
densityPlot(megdata)
mdsPlot(megdata, sampGroups =megdata@colData$Sex, sampNames = megdata@colData$Sample_ID)
megdata@colData$Sex<-ifelse(megdata@colData$Sample_ID=="TS222", "F", 
                            ifelse(megdata@colData$Sample_ID=="TS247", "M",megdata@colData$Sex))
mdsPlot(megdata, sampGroups =megdata@colData$Sex, sampNames = megdata@colData$Sample_ID)
mdsPlot(megdata, sampGroups =megdata@colData$Slide, sampNames = megdata@colData$Sample_ID)


#plot controls
plotCtrl(megdata,IDorder=NULL)

#New threshold 0.1 (10% of CpGs low call in sample)
qcscore2<-QCinfo(megdata, detPthre=0.000001, nbthre=3, CpGthre=0.05, samplethre=0.05,outlier=TRUE, distplot=T)
pData(epicData)$Bad4<-ifelse(rownames(pData(epicData))%in%qcscore2$badsample, 1,0)

phenomeg<-as.data.frame(pData(megdata))
levels(as.factor(cordreferencefull$CellType))

cell.types<-c("Gran","Bcell", "CD4T" , "CD8T" ,  "Mono" , "NK"  ,  "nRBC")

boxplot(phenomeg[,cell.types])



#IDOL optimized


#Step 1 IDOL
#load("D:/OneDrive/OneDrive - Dartmouth College/TReg IDOL/Treg_Reference_library.RData")
library(IDOL) #This functions will be intgrated in the FlowSorted.Blood.EPIC package in the near future


Reference2<-preprocessNoob(cordreferencefull)
mdsPlot(Reference2, sampGroups = Reference2$CellType)


cellTypes <- c("Gran","Bcell", "CD4T" , "CD8T" ,  "Mono" , "NK"  ,  "nRBC")






referenceBetas<-getBeta(Reference2)
referenceCovars<-as.data.frame(pData(Reference2))
candidates<-CandidateDMRFinder.v2(cellTypes, referenceBetas, referenceCovars, M = 1000,
                                  equal.variance = F)




betasmeg<-preprocessNoob(megdata)
allBetasmix<-getBeta(betasmeg)
phenomix<-phenomeg


phenotrain<-phenomix[1:12, ]
phenotest<-phenomix[13:24, ]
trainingBetas<-allBetasmix[,1:12]
testingBetas<-allBetasmix[,13:24]

set.seed(1234) 
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 150, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 200, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 250, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 300, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 350, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 400, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 450, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 500, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 550, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 600, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 650, maxIt = 500, numCores = 20)
IDOLoptimize(candidates, trainingBetas, phenotrain,
             libSize = 700, maxIt = 500, numCores = 20)
#finalIDOL750<-IDOLoptimize(candidates, trainingBetas, phenotrain,
#                           libSize = 750, maxIt = 500, numCores = 20)
#system.time(finalIDOL800<-IDOLoptimize(candidates, trainingBetas, phenotrain,
#                                       libSize = 800, maxIt = 500, numCores = 20))

#Compare the libraries testing vs training 
librarycomp<-matrix(data = NA, nrow = 12, ncol = 4, 
                    dimnames = list(c("150", "200", "250","300","350", "400", "450", "500", "550", "600", "650", "700"),
                                    c("R2_training", "RMSE_training", "R2_testing", "RMSE_testing")))

R2compute = function(obs, pred) {
    r2 = NULL
    for(i in 1:dim(obs)[2]) {
        y = obs[,i]
        x = pred[,i]
        r2[i] = summary(lm(y~x))[[8]]
    }
    sum(r2)/dim(obs)[2]
}

librarylist<-c("150", "200", "250", "300","350", "400", "450", "500", "550", "600", "650", "700")

for (i in librarylist){
    load(file = paste("IDOL optimized DMR library_", i, ".RData", sep = ""))
    ctpred = minfi:::projectCellType(trainingBetas[IDOL.optim.DMRs,], IDOL.optim.coefEsts)
    omega.tilde = 100*ctpred
    omega.obs = phenotrain[,cellTypes]
    # compute RMSE based on all probes
    # Whole blood
    diff = as.matrix(omega.tilde - omega.obs)
    RMSE_training = sqrt(sum(diff^2)/length(cellTypes))
    # compute R2 based on all probes
    R2_training = R2compute(omega.obs, omega.tilde)
    ctpred = minfi:::projectCellType(testingBetas[IDOL.optim.DMRs,], IDOL.optim.coefEsts)
    omega.tilde = 100*ctpred
    omega.obs = phenotest[,cellTypes]
    # compute RMSE based on all probes
    # Whole blood
    diff = as.matrix(omega.tilde - omega.obs)
    RMSE_testing = sqrt(sum(diff^2)/length(cellTypes))
    # compute R2 based on all probes
    R2_testing = R2compute(omega.obs, omega.tilde)
    librarycomp[i,]<-c(R2_training, RMSE_training, R2_testing,RMSE_testing)
}

librarycompcord<-librarycomp

library(reshape2)
librarycompR2<-melt(librarycomp[,c(1,3)])
librarycompRMSE<-melt(librarycomp[,c(2,4)])
library(ggplot2)
ggplot(data=librarycompR2, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity", position=position_dodge())+labs(x="Library size", y=expression(R^2), fill=NULL)#+  scale_y_continuous(limits =c(0.90, 1))
ggplot(data=librarycompRMSE, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity", position=position_dodge())+labs(x="Library size", y="RMSE", fill=NULL)


#550 was better in any testing and training
load("D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive/IDOL cord/IDOL optimized DMR library_550.RData")
length(IDOL.optim.DMRs)
library(pheatmap)
pheatmap(IDOL.optim.coefEsts,color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, annotation_col = annotation_col, 
         fontsize = 14,
         show_rownames = F, show_colnames = T, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = paste(length(IDOL.optim.DMRs),"probes IDOL library"))

# Function that returns Root Mean Squared Error
rmse <- function(error)
{
    sqrt(mean(error^2))
}
coefs<-IDOL.optim.coefEsts
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno<-as.data.frame(anno)
coefs<-coefs[order(rownames(coefs)),]
anno<-anno[order(rownames(anno)),]
anno2<-anno[rownames(anno)%in% rownames(coefs),]
identical(rownames(anno2), rownames(coefs))
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno<-as.data.frame(anno)
coefs<-coefs[order(rownames(coefs)),]
anno<-anno[order(rownames(anno)),]
anno2<-anno[rownames(anno)%in% rownames(coefs),]
identical(rownames(anno2), rownames(coefs))
anno2$DNase<-ifelse(anno2$DHS==TRUE, "Yes", "No")
anno2$Enhancer<-ifelse(anno2$Enhancer==TRUE, "Yes", "No")
anno2$chrX<-ifelse(anno2$chr=="chrX", "Yes", "No")
annotation_row<-data.frame(DHS=anno2$DNase, Enhancer=anno2$Enhancer, ChrX=anno2$chrX, row.names = rownames(anno2))
annotation_col<-data.frame(CellType=colnames(coefs), row.names = colnames(coefs))
table(anno2$DNase)
table(anno2$Enhancer)
#table(anno2$Methyl450_Loci)#142 in EPIC array
table(anno2$Relation_to_Island)
anno4<-strsplit(anno2$UCSC_RefGene_Group, ";")
#character 0 should be transformed to NA or the lapply will fail
#anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) NA_character_ else x)
anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) "Intergenic" else x)
anno5<-lapply(anno5a, `[[`, 1)
anno6<-factor(unlist(anno5), levels= c("TSS1500", "TSS200", "1stExon", "5'UTR","Body", "3'UTR", "Intergenic"))
table(anno6)


anno7<-strsplit(anno2$UCSC_RefGene_Name, ";")
#character 0 should be transformed to NA or the lapply will fail
#anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) NA_character_ else x)
anno8a<-lapply(anno7, function(x) if(identical(x, character(0))) "Intergenic" else x)
anno8<-lapply(anno8a, `[[`, 1)
anno9<-as.character(unlist(anno8))
table(anno9)
anno2$gene<-anno9

ann_colors = list(CellType=c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3", Mono="#E7298A",
                             Gran="#66A61E", NK="#E6AB02", nRBC="red"),
                  Enhancer=c(No="black", Yes="lightgray"))



pheatmap(coefs,color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         annotation_row= annotation_row, annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = T, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "550 probes IDOL selection cord blood")


library(FDb.InfiniumMethylation.hg19)
#hm450 <- get450k()
load("D:/OneDrive/OneDrive - Dartmouth College/Crossreactive/hm450.manifest.rda")
neartss<- getNearestTSS(hm450.manifest)
neargene<-getNearestGene(hm450.manifest)
neartss2<-neartss[rownames(neartss)%in%rownames(anno2),]
neargene2<-neargene[rownames(neargene)%in%rownames(anno2),]
anno2$function_group<-anno6
annomask<-as.data.frame(hm450.manifest)
annomask<-annomask[rownames(annomask)%in% rownames(anno2),]

general_anno<-cbind.data.frame(coefs, anno2[order(rownames(anno2)),], neartss2[order(rownames(neartss2)),],neargene2[order(rownames(neargene2)),],
                               annomask[order(rownames(annomask)),])
write.csv(general_anno, file="IDOL CpGs cord ref only 550.csv")



celltypes<-c( "CD4T","CD8T","Bcell", "NK","Gran","Mono", "nRBC")

truevalues<-phenomix[,celltypes]
counts1 = minfi:::projectCellType(allBetasmix[IDOL.optim.DMRs,], IDOL.optim.coefEsts)
#counts1<-countsEPICTreg$counts#Use the estimates below

celltypes2<-rev(celltypes)
counts1<-as.data.frame(counts1*100)
counts1<-counts1[,celltypes]
difcounts1<-counts1-truevalues

difcounts1<-difcounts1[,rev(colnames(difcounts1))]

colorcelltype2<-c(nRBC="red", Mono="#E7298A",
                  Gran="#66A61E", NK="#E6AB02", Bcell="#1B9E77",CD8T ="#7570B3", CD4T="#D95F02")

colorcelltype<-rev(colorcelltype2)

boxplot(difcounts1,  horizontal=T, col=colorcelltype2)#, yaxt = "n", add = TRUE, boxwex=0.25, border=c("blue", "blue","blue","blue","blue","blue"),#boxfill="red",
#at = 1:ncol(difall[,-7]) - 0.225) #shift these left by -0.15





#Cord blood only 
pheno2<-phenomix
truevalues<-phenomix[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(25:35) ]<-"Training"
pheno2$clas[c(36:37,39:43, 45:48) ]<-"Testing"
#pheno2$clas[c(25:36) ]<-"Training"#, 25:36) ]<-"Training"
#pheno2$clas[c(37:48) ]<-"Testing"#, 37:48) ]<-"Testing"
pheno2<-pheno2#[-14,]
cellcount_EPICIDOL<-counts1#[-c(14,20),]
truevalues<-truevalues#[-c(14,20),]
rownames(cellcount_EPICIDOL)<-rownames(truevalues)


par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[c(25:37, 39:43, 45:48),i]~truevalues[c(25:37, 39:43, 45:48),i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[c(25:37, 39:43, 45:48),i]~truevalues[c(25:37, 39:43, 45:48),i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[c(25:37, 39:43, 45:48),i]~truevalues[c(25:37, 39:43, 45:48),i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}


#Cord blood only training
pheno2<-phenomix
truevalues<-phenomix[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(25:36) ]<-"Training"#, 25:36) ]<-"Training"
pheno2$clas[c(37:48) ]<-"Testing"#, 37:48) ]<-"Testing"
pheno2<-pheno2#[-14,]
cellcount_EPICIDOL<-counts1#[-14,]
cellcount_EPICIDOL<-round(cellcount_EPICIDOL,1)

truevalues<-truevalues#[-14,]
rownames(cellcount_EPICIDOL)<-rownames(truevalues)
par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        #par(new=TRUE)
        #plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
        #     cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
        #     xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
        #     xlim=c(0,80), ylim=c(0,80),
        #     pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
        #     cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[25:36,i]~truevalues[25:36,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        #par(new=TRUE)
        #plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
         #    cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
          #   xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
           #  xlim=c(0,50), ylim=c(0,50),
            # pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             #cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[25:36,i]~truevalues[25:36,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        #par(new=TRUE)
        #plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
         #    cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
          #   xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
           #  xlim=c(0,30), ylim=c(0,30),
            # pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             #cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[25:36,i]~truevalues[25:36,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}

#Cord blood only Testing
pheno2<-phenomix
truevalues<-phenomix[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(25:36) ]<-"Training"#, 25:36) ]<-"Training"
pheno2$clas[c(37:48) ]<-"Testing"#, 37:48) ]<-"Testing"
pheno2<-pheno2#[-14,]
cellcount_EPICIDOL<-counts1#[-14,]
cellcount_EPICIDOL<-round(cellcount_EPICIDOL,1)

truevalues<-truevalues#[-14,]
rownames(cellcount_EPICIDOL)<-rownames(truevalues)
par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        #plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
         #    cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
          #   xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
           #  xlim=c(0,80), ylim=c(0,80),
            # pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             #cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        #par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[37:48,i]~truevalues[37:48,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        #plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
         #    cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
          #   xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
           #  xlim=c(0,50), ylim=c(0,50),
            # pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             #cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        #par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[37:48,i]~truevalues[37:48,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        #plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
         #    cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
          #   xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
           #  xlim=c(0,30), ylim=c(0,30),
            # pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             #cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        #par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[37:48,i]~truevalues[37:48,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}


#pheatmap

allrefnoob2<-allrefnoob[, !(allrefnoob$CellType=="MIX" & allrefnoob$Study!="Salas")]

anno7<-strsplit(anno2$UCSC_RefGene_Name, ";")
#character 0 should be transformed to NA or the lapply will fail
#anno5a<-lapply(anno4, function(x) if(identical(x, character(0))) NA_character_ else x)
anno8a<-lapply(anno7, function(x) if(identical(x, character(0))) "Intergenic" else x)
anno8<-lapply(anno8a, `[[`, 1)
anno9<-as.character(unlist(anno8))
table(anno9)
anno2$gene<-anno9
anno2$CD8<-ifelse(anno2$gene=="CD8A" |anno2$gene=="CD8B", "Yes", "No")
annotation_row<-data.frame(DHS=anno2$DNase, Enhancer=anno2$Enhancer, ChrX=anno2$chrX, CD8=anno2$CD8, row.names = rownames(anno2))



annotation_col<-data.frame(CellType=allrefnoob2$CellType, CellType_old=allrefnoob2$CellType2, Study= allrefnoob2$Study, row.names = colnames(allrefnoob2))

ann_colors = list(CellType=c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3", Eos= "pink", Gran="#66A61E",
                             MIX="black", Mono="#E7298A",Neu="lightgreen",  NK="#E6AB02", nRBC="red", 
                             PanT="darkblue", PBMC="lightgray",  WBC="darkgray"),
                  CellType_old=c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3", Eos= "pink", Gran="#66A61E",
                             MIX="black", Mono="#E7298A",Neu="lightgreen",  NK="#E6AB02", nRBC="red", 
                             PanT="darkblue", PBMC="lightgray",  WBC="darkgray"),
                  Study=c( "Bakulski"="blue", "Gervin"="red", "deGoede"="black", Reinius="yellow", Salas= "Green"),
                  Enhancer=c(No="black", Yes="lightgray"),  CD8=c(No="black", Yes="lightgray"))


betascells<-getBeta(allrefnoob2)
betascells<-betascells[rownames(betascells)%in% IDOL.optim.DMRs,]

pheatmap(betascells,color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         annotation_row= annotation_row, annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "550 probes IDOL selection adult+cord blood combined referenxw")


pheatmap(betascells[rownames(betascells)%in%rownames(anno2[anno2$CD8=="Yes",]),],color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, 
         annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "550 probes IDOL selection adult+cord blood combined referenxw")



library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]
FlowSorted.Blood.EPIC

allrefnew#use make-data.R
allrefnew$CellType2<-allrefnew$CellType_original
allrefnew$Cord<-"Yes"

allrefnewnoob<-preprocessNoob(allrefnew)

mdsPlot(allrefnewnoob[, allrefnewnoob$CellType!="Endothelial" & allrefnewnoob$CellType!="Epithelial" & allrefnewnoob$CellType!="Stromal"], 
        sampGroups = allrefnewnoob$CellType[allrefnewnoob$CellType!="Endothelial" & allrefnewnoob$CellType!="Epithelial" & allrefnewnoob$CellType!="Stromal"], 
        pal = c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3",  Eos="pink",#Endothelial="darkblue", Epithelial= "maroon", 
                Gran="#66A61E", MIX="black", Mono="#E7298A", Neu= "lightgreen",
                NK="#E6AB02", nRBC="red", PanT="darkblue", PBMC="darkgray", WBC="gray" ), 
        main= "Top 1000 most variable probes restricted to blood cells")

mdsPlot(allrefnewnoob[IDOLOptimizedCpGs450klegacy, allrefnewnoob$CellType!="Endothelial" & allrefnewnoob$CellType!="Epithelial" & allrefnewnoob$CellType!="Stromal"], 
        sampGroups = allrefnewnoob$CellType[allrefnewnoob$CellType!="Endothelial" & allrefnewnoob$CellType!="Epithelial" & allrefnewnoob$CellType!="Stromal"], 
        pal = c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3",  Eos="pink",#Endothelial="darkblue", Epithelial= "maroon", 
                Gran="#66A61E", MIX="black", Mono="#E7298A", Neu= "lightgreen",
                NK="#E6AB02", nRBC="red", PanT="darkblue", PBMC="darkgray", WBC="gray" ), 
        main= "350 IDOL probes restricted to blood cells")

allrefnewnoob$CellType2<-ifelse(is.na(allrefnewnoob$CellType2), allrefnewnoob$CellType, allrefnewnoob$CellType2)
allrefnewnoob$Cord<-ifelse(allrefnewnoob$Study=="Reinius" | allrefnewnoob$Study=="Salas", "No", "Yes")

allrefnewnoob2<-allrefnewnoob[, allrefnewnoob$CellType!="Endothelial" & allrefnewnoob$CellType!="Epithelial" & allrefnewnoob$CellType!="Stromal"]

save(allrefnew, allrefnewnoob, file="All 450 EPIC references.RData")

ann_colors = list(CellType=c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3",  Eos="pink",#Endothelial="darkblue", Epithelial= "maroon", 
                             Gran="#66A61E", MIX="black", Mono="#E7298A", Neu= "lightgreen",
                             NK="#E6AB02", nRBC="red", PanT="darkblue", PBMC="darkgray", WBC="gray" ),
                  CellType2=c(Bcell="#1B9E77", CD4T="#D95F02", CD8T ="#7570B3",  Eos="pink",#Endothelial="darkblue", Epithelial= "maroon", 
                             Gran="#66A61E", MIX="black", Mono="#E7298A", Neu= "lightgreen",
                             NK="#E6AB02", nRBC="red", PanT="darkblue", PBMC="darkgray", WBC="gray" ),
                  Enhancer=c(No="black", Yes="lightgray"))

annotation_col=data.frame(CellType=allrefnewnoob2$CellType, CellType2=allrefnewnoob2$CellType2, Cord= allrefnewnoob2$Cord, Study= allrefnewnoob2$Study,#Ethnic=as.factor(noobcordEPIC$ethnicity), 
                          row.names = colnames(allrefnewnoob2))


pheatmap(getBeta(allrefnewnoob2)[IDOLOptimizedCpGs450klegacy,],color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, 
         annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "350 probes IDOL adult selection adult cord blood")

load("D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive/IDOL cord/IDOL optimized DMR library_550.RData")

pheatmap(getBeta(allrefnewnoob2)[rownames(allrefnewnoob2)%in%IDOL.optim.DMRs,],color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, 
         annotation_col = annotation_col, fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "517/550 probes IDOL cord blood selection, adult cord blood")


pheatmap(getBeta(allrefnewnoob2)[rownames(allrefnewnoob2)%in%IDOL.optim.DMRs, allrefnewnoob2$Cord=="Yes"],color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, 
         annotation_col = annotation_col[allrefnewnoob2$Cord=="Yes",], fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "517/550 probes IDOL cord blood selection, cord blood only")


pheatmap(getBeta(allrefnewnoob2)[rownames(allrefnewnoob2)%in%IDOL.optim.DMRs, allrefnewnoob2$Cord=="Yes" & allrefnewnoob2$CellType!="MIX"],
         color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, 
         annotation_col = annotation_col[allrefnewnoob2$Cord=="Yes" & allrefnewnoob2$CellType!="MIX",], fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "517/550 probes IDOL cord blood selection, cord blood only clean")

pheatmap(getBeta(allrefnewnoob2)[IDOLOptimizedCpGs450klegacy, allrefnewnoob2$Cord=="Yes" & allrefnewnoob2$CellType!="MIX"],
         color = colorRampPalette(c("yellow", "black", "blue"), space = "Lab")(128),
         #annotation_row= annotation_row, 
         annotation_col = annotation_col[allrefnewnoob2$Cord=="Yes" & allrefnewnoob2$CellType!="MIX",], fontsize = 14,
         show_rownames = F, show_colnames = F, legend = T, cluster_rows = T, annotation_colors = ann_colors,
         main = "350 probes IDOL adult blood selection, cord blood only clean")


allcordclean<-allrefnew[, allrefnew$Cord=="Yes" & allrefnew$CellType!="MIX" & allrefnew$CellType!="Endothelial" & allrefnew$CellType!="Epithelial" & allrefnew$CellType!="Stromal"]
table(allcordclean$CellType)

estimatescordrefadultIDOL<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                       processMethod = "preprocessNoob",
                                       probeSelect = "IDOL", 
                                       cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                     "Mono", "Gran", "nRBC"), 
                                       referencePlatform = 
                                           "IlluminaHumanMethylationEPIC",
                                       referenceset = "allcordclean",
                                       IDOLOptimizedCpGs =IDOLOptimizedCpGs450klegacy, 
                                       returnAll = F)

pheno2<-phenomeg
truevalues<-phenomeg[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(1:24) ]<-"Training"
cellcount_EPICIDOL<-estimatescordrefadultIDOL$counts*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

pheno2<-pheno2[-c(14,20),]
cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
truevalues<-truevalues[-c(14,20),]

par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}



estimatescordrefcordIDOL<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                               processMethod = "preprocessNoob",
                                               probeSelect = "IDOL", 
                                               cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                             "Mono", "Gran", "nRBC"), 
                                               referencePlatform = 
                                                   "IlluminaHumanMethylationEPIC",
                                               referenceset = "allcordclean",
                                               IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                               returnAll = F)

pheno2<-phenomeg
truevalues<-phenomeg[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(1:24) ]<-"Training"
cellcount_EPICIDOL<-estimatescordrefcordIDOL$counts*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

pheno2<-pheno2[-c(14,20),]
cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
truevalues<-truevalues[-c(14,20),]


par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}


#Leave one out
#Lin

Study<-levels(as.factor(allcordclean$Study))

for (j in Study){
    allcordclean1<-allcordclean[,allcordclean$Study!=j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "IDOL", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran", "nRBC"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}


#####
#Leave one out as originally published
####


allcordunclean<-allrefnew[, allrefnew$Cord=="Yes" & allrefnew$CellType!="Endothelial" & allrefnew$CellType!="Epithelial" & allrefnew$CellType!="Stromal"]
table(allcordunclean$CellType, allcordunclean$CellType2)

allcordunclean$CellType3<-allcordunclean$CellType
allcordunclean$CellType<-ifelse(allcordunclean$CellType2=="nRBC" & allcordunclean$CellType3=="MIX", "MIX", allcordunclean$CellType2)
table(allcordunclean$CellType, allcordunclean$CellType2)

table(allcordunclean$CellType, allcordunclean$Study)


estimatescordrefcordIDOL<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                              processMethod = "preprocessNoob",
                                              probeSelect = "IDOL", 
                                              cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                            "Mono", "Gran", "nRBC"), 
                                              referencePlatform = 
                                                  "IlluminaHumanMethylationEPIC",
                                              referenceset = "allcordunclean",
                                              IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                              returnAll = F)

pheno2<-phenomeg
truevalues<-phenomeg[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(1:24) ]<-"Training"
cellcount_EPICIDOL<-estimatescordrefcordIDOL$counts*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

pheno2<-pheno2[-c(14,20),]
cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
truevalues<-truevalues[-c(14,20),]


par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}


#Leave one out
#Lin

Study<-levels(as.factor(allcordunclean$Study))

for (j in Study){
    allcordclean1<-allcordunclean[,allcordunclean$Study!=j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "IDOL", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran", "nRBC"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}



#One by one uncleaned


Study<-levels(as.factor(allcordunclean$Study))

for (j in Study[1:2]){
    allcordclean1<-allcordunclean[,allcordunclean$Study==j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "IDOL", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran", "nRBC"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}



for (j in Study[3:4]){
    allcordclean1<-allcordunclean[,allcordunclean$Study==j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "IDOL", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes[1:6]) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}

#One by one cleaned


Study<-levels(as.factor(allcordclean$Study))

for (j in Study[1:2]){
    allcordclean1<-allcordclean[,allcordclean$Study==j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "IDOL", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran", "nRBC"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}



for (j in Study[3:4]){
    allcordclean1<-allcordclean[,allcordclean$Study==j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "IDOL", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =IDOL.optim.DMRs, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes[1:6]) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}




#One by one minfi  uncleaned


Study<-levels(as.factor(allcordunclean$Study))

for (j in Study[1:2]){
    allcordclean1<-allcordunclean[,allcordunclean$Study==j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "any", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran", "nRBC"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =NULL, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}



for (j in Study[3:4]){
    allcordclean1<-allcordunclean[,allcordunclean$Study==j]
    
    estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                                   processMethod = "preprocessNoob",
                                                   probeSelect = "any", 
                                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                                 "Mono", "Gran"), 
                                                   referencePlatform = 
                                                       "IlluminaHumanMethylationEPIC",
                                                   referenceset = "allcordclean1",
                                                   IDOLOptimizedCpGs =NULL, 
                                                   returnAll = F)
    
    pheno2<-phenomeg
    truevalues<-phenomeg[,celltypes]
    pheno2$clas<-NA
    
    pheno2$clas[c(1:24) ]<-"Training"
    cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    pheno2<-pheno2[-c(14,20),]
    cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
    truevalues<-truevalues[-c(14,20),]
    
    
    
    par(mfrow=c(2,4))
    for (i in celltypes[1:6]) {
        #error <- truevalues[,i] - counts[,i]*100
        
        if (i=="Gran"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,80), ylim=c(0,80),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        } else if (i=="CD4T"){
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,50), ylim=c(0,50),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
        else {
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            par(new=TRUE)
            plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
                 cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
                 xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
                 xlim=c(0,30), ylim=c(0,30),
                 pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
                 cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
            lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
            rmse(lm.D9$residuals)
            summary(lm.D9)$r.squared
            abline(lm.D9, col=colorcelltype[i])
            abline(a = 0, b=1, col="black", lwd=1, lty=2)
            text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
            text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
            
        }
    }   
}

#All ref minfi

allcordclean1<-allcordunclean

estimatescordrefcordIDOL1<-estimateCellCounts2(megdata, compositeCellType = "Blood", 
                                               processMethod = "preprocessNoob",
                                               probeSelect = "any", 
                                               cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                             "Mono", "Gran", "nRBC"), 
                                               referencePlatform = 
                                                   "IlluminaHumanMethylationEPIC",
                                               referenceset = "allcordclean1",
                                               IDOLOptimizedCpGs =NULL, 
                                               returnAll = F)

pheno2<-phenomeg
truevalues<-phenomeg[,celltypes]
pheno2$clas<-NA

pheno2$clas[c(1:24) ]<-"Training"
cellcount_EPICIDOL<-estimatescordrefcordIDOL1$counts*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

pheno2<-pheno2[-c(14,20),]
cellcount_EPICIDOL<-cellcount_EPICIDOL[-c(14,20),]
truevalues<-truevalues[-c(14,20),]



par(mfrow=c(2,4))
for (i in celltypes) {
    #error <- truevalues[,i] - counts[,i]*100
    
    if (i=="Gran"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i, "no", j),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i, "no", j),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,80), ylim=c(0,80),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(15, 65, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(15, 55, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    } else if (i=="CD4T"){
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,50), ylim=c(0,50),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(10, 45, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(10, 37, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
    else {
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Training",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Training",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 21,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        par(new=TRUE)
        plot(truevalues[rownames(truevalues)%in% rownames(pheno2[pheno2$clas=="Testing",]),i],
             cellcount_EPICIDOL[rownames(cellcount_EPICIDOL)%in% rownames(pheno2[pheno2$clas=="Testing",]),i], main=paste(i),
             xlab="",ylab="",#paste("True", i, "%"), ylab=paste("Estimated", i, "%"),
             xlim=c(0,30), ylim=c(0,30),
             pch = 22,  col="black", bg=colorcelltype[i], #pch = 20,#as.numeric(mix),
             cex=2.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        abline(lm.D9, col=colorcelltype[i])
        abline(a = 0, b=1, col="black", lwd=1, lty=2)
        text(5, 27, substitute(paste(R^2,m), list(m=round(summary(lm.D9)$r.squared[1]*100,1))), cex = 1.5)
        text(5, 23, substitute(paste("RMSE",m),list(m=round(rmse(lm.D9$residuals)[1],2))), cex = 1.5)
        
    }
}   




Reference4<-preprocessNoob(allcordclean[,allcordclean$CellType%in%cellTypes])

cellIdentity = Reference4$CellType  # identity of the cell types 

# load in the IDOL optimized library (size = 400 CpGs) for the deconvoluting mixtures consisting of the normal 6 
# leukocyte subtypes
#load("IDOL Optimization_400.RData")
betas.normal6<-getBeta(Reference4)


# extract the IDOL optimized CpGs
IDOLOptim517 = IDOL.optim.DMRs[IDOL.optim.DMRs%in%rownames(Reference4)]#IDOLObjects[["IDOL Optimized Library"]]


#minfi process
#Original distribution
Reference5<-allcordunclean[, allcordunclean$Cord=="Yes" & allcordunclean$CellType!="PanT" & allcordunclean$CellType!="MIX" & allcordunclean$CellType!="WBC"]


#Cell contamination
phenounclean<-as.data.frame(pData(allcordunclean))

countsunclean<-estimateCellCounts2(allcordunclean, compositeCellType = "Blood", 
                                   processMethod = "preprocessNoob",
                                   probeSelect = "IDOL", 
                                   cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                                 "Mono", "Neu"), 
                                   referencePlatform = 
                                       "IlluminaHumanMethylationEPIC",
                                   referenceset = NULL,
                                   IDOLOptimizedCpGs =IDOLOptimizedCpGs450klegacy, 
                                   returnAll = FALSE)

head(countsunclean$counts)
pheno3<-round(countsunclean$counts*100,2)
pheno3<-as.data.frame(pheno3)
pheno3$total<-rowSums(pheno3)
pheno3a<-pheno3[,c("CD8T", "CD4T", "NK", "Bcell", 
                   "Mono", "Neu")]*100/pheno3$total
pheno3a$total<-rowSums(pheno3a)#To confirm it is 100%
pheno3<-pheno3a[,c("CD8T", "CD4T", "NK", "Bcell", 
                   "Mono", "Neu")]

pheno3$CellType<-phenounclean$CellType
pheno3$CellType2<-phenounclean$CellType3
pheno3$Study<-phenounclean$Study
pheno3$sampname<-paste0(phenounclean$CellType,"_",phenounclean$Study,"_", seq(1:dim(pheno3)[1]) )
head(pheno3)
pheno3<-pheno3[order(pheno3$sampname),]
pheno3$sampnameorig<-rownames(pheno3)
rownames(pheno3)<-pheno3$sampname

library(ggplot2)
library(reshape)
attach(pheno3)
stack<-pheno3[,c("CD8T", "CD4T", "NK", "Bcell", 
                                         "Mono", "Neu")]
stack<-t(stack)
dat2 <- melt(stack)#, id.vars = "Sample")

fetaladultspine<- ggplot(dat2, aes(x=X2, y=value, fill=X1)) + 
    geom_bar(stat="identity") +
    xlab("\nSample") +
    ylab("Cell %\n") +
    #guides(fill=T) +
    theme_bw()+
    ggtitle("Cell projection IDOL") + scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3",  
                                                                 "#E7298A", "#66A61E",
                                                                 "#E6AB02", "red"))+
    theme(axis.text.x=element_text(angle = -90, hjust = 0))+
    theme(legend.title=element_blank()) 
#theme(legend.position="none")
#detach(pheno2)
fetaladultspine
pheno3$CellType<-ifelse(pheno3$CellType=="Neu", "Gran", pheno3$CellType)
pheno3$Gran=pheno3$Neu
#Bcell
library(forcats)
for (i in c("CD8T", "CD4T", "NK", "Bcell","Mono", "Gran")){#, "nRBC")){
    stack<-pheno3[pheno3$CellType==i,c("CD8T", "CD4T", "NK", "Bcell", 
                                             "Mono", "Gran")]
    stack<-t(stack)
    dat2 <- melt(stack)#, id.vars = "Sample")
    dat2$X1 <- relevel(dat2$X1, ref=i)
    levels(dat2$X1)
    dat2$X1 <-fct_rev(dat2$X1)
    levels(dat2$X1)
    colors1<-c(Bcell= "#1B9E77", CD4T="#D95F02", CD8T= "#7570B3",  
               Mono= "#E7298A", Gran= "#66A61E",
               NK= "#E6AB02")#, nRBC="red")
    colors2<-colors1[levels(dat2$X1)]
    fetaladultspine<- ggplot(dat2, aes(x=X2, y=value, fill=X1)) + 
        geom_bar(stat="identity") +
        xlab("") +#xlab("\nSample") +
        ylab("") +#ylab("Cell %\n") + 
        geom_hline(yintercept=70, linetype="dashed", color = "red")+ 
        #guides(fill=T) +
        theme_bw()+
        ggtitle(paste("Cell projection IDOL", i)) + scale_fill_manual(values=colors2)+
        theme(axis.text.x=element_text(angle = -90, hjust = 0))+
        theme(legend.title=element_blank(), legend.position = "none") 
    #theme(legend.position="none")
    #detach(pheno2)
    assign(x = paste0(i,"_bar"),value = fetaladultspine)   
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}


multiplot(  Mono_bar, NK_bar, CD8T_bar, CD4T_bar, Gran_bar,Bcell_bar, cols=3)



allcordcleanstr<-allcordclean[,-c(which(allcordclean$CellType=="CD8T" & allcordclean$CellType2=="CD4T"),
                                 which(allcordclean$CellType=="NK" & allcordclean$CellType2=="CD8T"),
                                 which(allcordclean$CellType=="CD4T" & allcordclean$CellType2=="CD8T"))]
table(allcordcleanstr$Study, allcordcleanstr$CellType)
table(allcordcleanstr$CellType, allcordcleanstr$CellType2)



#Quick loop
library(FlowSorted.Blood.EPIC)
Study<-levels(as.factor(allcordclean$Study))
allcordcleannoob<-preprocessNoob(allcordclean)
megdatanoob<-preprocessNoob(megdata)
megdatanoob<-megdatanoob[rownames(megdatanoob)%in%rownames(allcordcleannoob),]
truevalues<-phenomeg[-c(14,20),celltypes]
write.csv(truevalues, file="True_values.csv")


Bakulski<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="Bakulski"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Bakulski$coefEsts)
deGoede<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="deGoede"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(deGoede$coefEsts)
Gervin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="Gervin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Gervin$coefEsts)
Lin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="Lin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Lin$coefEsts)
Allref<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Allref$coefEsts)

cleancoefs<-list(Bakulski=Bakulski$coefEsts, deGoede=deGoede$coefEsts, Gervin=Gervin$coefEsts, Lin=Lin$coefEsts, Allref=Allref$coefEsts)



cleanminfi<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                "deGoede_R2", "deGoede_rmse",
                                                                                "Gervin_R2", "Gervin_rmse",
                                                                                "Lin_R2", "Lin_rmse",
                                                                                "Allref_R2", "Allref_rmse")))


for (j in Study[1:2]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(cleancoefs[[j]]),-c(14,20)], coefCellType = cleancoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("minfi_cleaned_",j,".csv"))
    for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    cleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    cleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}
}

for (j in Study[3:4]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(cleancoefs[[j]]),-c(14,20)], coefCellType = cleancoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("minfi_cleaned_",j,".csv"))
    for(i in celltypes[1:6]){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
truevalues<-phenomeg[-c(14,20),celltypes]

estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(cleancoefs[[j]]),-c(14,20)], coefCellType = cleancoefs[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    cleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    cleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}


#Unclean
#Quick loop
library(FlowSorted.Blood.EPIC)
Study<-levels(as.factor(allcordclean$Study))
allcorduncleannoob<-preprocessNoob(allcordunclean)
megdatanoob<-preprocessNoob(megdata)
megdatanoob<-megdatanoob[rownames(megdatanoob)%in%rownames(allcorduncleannoob),]

Bakulski<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcorduncleannoob[,allcorduncleannoob$Study=="Bakulski"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Bakulski$coefEsts)
deGoede<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcorduncleannoob[,allcorduncleannoob$Study=="deGoede"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(deGoede$coefEsts)
Gervin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcorduncleannoob[,allcorduncleannoob$Study=="Gervin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Gervin$coefEsts)
Lin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcorduncleannoob[,allcorduncleannoob$Study=="Lin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Lin$coefEsts)
Allref<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcorduncleannoob, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Allref$coefEsts)

uncleancoefs<-list(Bakulski=Bakulski$coefEsts, deGoede=deGoede$coefEsts, Gervin=Gervin$coefEsts, Lin=Lin$coefEsts, Allref=Allref$coefEsts)

uncleanminfi<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                "deGoede_R2", "deGoede_rmse",
                                                                                "Gervin_R2", "Gervin_rmse",
                                                                                "Lin_R2", "Lin_rmse",
                                                                                "Allref_R2", "Allref_rmse")))


for (j in Study[1:2]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(uncleancoefs[[j]]),-c(14,20)], coefCellType = uncleancoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    for(i in celltypes){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        uncleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        uncleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

for (j in Study[3:4]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(uncleancoefs[[j]]),-c(14,20)], coefCellType = uncleancoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    
    for(i in celltypes[1:6]){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        uncleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        uncleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
truevalues<-phenomeg[-c(14,20),celltypes]

estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(uncleancoefs[[j]]),-c(14,20)], coefCellType = uncleancoefs[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    uncleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    uncleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}


#IDOL
load("D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive/IDOL cord/IDOL optimized DMR library_550.RData")
IDOLOptim517<-IDOL.optim.DMRs[IDOL.optim.DMRs%in%rownames(megdatanoob)]
Bakulski<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="Bakulski"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Bakulski$compTable[IDOLOptim517,3:9])
deGoede<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="deGoede"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(deGoede$compTable[IDOLOptim517,3:9])
Gervin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="Gervin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Gervin$compTable[IDOLOptim517,3:8])
Lin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study=="Lin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Lin$compTable[IDOLOptim517,3:8])
Allref<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Allref$compTable[IDOLOptim517,3:9])

IDOLcoefs<-list(Bakulski=as.matrix(Bakulski$compTable[IDOLOptim517,3:9]), deGoede=as.matrix(deGoede$compTable[IDOLOptim517,3:9]), Gervin=as.matrix(Gervin$compTable[IDOLOptim517,3:8]), 
                Lin=as.matrix(Lin$compTable[IDOLOptim517,3:8]), Allref=as.matrix(Allref$compTable[IDOLOptim517,3:9]))

cleanIDOL<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                "deGoede_R2", "deGoede_rmse",
                                                                                "Gervin_R2", "Gervin_rmse",
                                                                                "Lin_R2", "Lin_rmse",
                                                                                "Allref_R2", "Allref_rmse")))


for (j in Study[1:2]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs[[j]]),-c(14,20)], coefCellType = IDOLcoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_",j,".csv"))
    for(i in celltypes){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

for (j in Study[3:4]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs[[j]]),-c(14,20)], coefCellType = IDOLcoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_",j,".csv"))
    for(i in celltypes[1:6]){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
truevalues<-phenomeg[-c(14,20),celltypes]

estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs[[j]]),-c(14,20)], coefCellType = IDOLcoefs[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    cleanIDOL[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    cleanIDOL[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}


#IDOL leave one out (loo)
Bakulski<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study!="Bakulski"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Bakulski$compTable[IDOLOptim517,3:9])
deGoede<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study!="deGoede"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(deGoede$compTable[IDOLOptim517,3:9])
Gervin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study!="Gervin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Gervin$compTable[IDOLOptim517,3:8])
Lin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob[,allcordcleannoob$Study!="Lin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Lin$compTable[IDOLOptim517,3:8])
Allref<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleannoob, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Allref$compTable[IDOLOptim517,3:9])

IDOLcoefs_loo<-list(Bakulski=as.matrix(Bakulski$compTable[IDOLOptim517,3:9]), deGoede=as.matrix(deGoede$compTable[IDOLOptim517,3:9]), Gervin=as.matrix(Gervin$compTable[IDOLOptim517,3:9]), 
                Lin=as.matrix(Lin$compTable[IDOLOptim517,3:9]), Allref=as.matrix(Allref$compTable[IDOLOptim517,3:9]))

cleanIDOL_loo<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                               "deGoede_R2", "deGoede_rmse",
                                                                               "Gervin_R2", "Gervin_rmse",
                                                                               "Lin_R2", "Lin_rmse",
                                                                               "Allref_R2", "Allref_rmse")))


for (j in Study){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs_loo[[j]]),-c(14,20)], coefCellType = IDOLcoefs_loo[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_",j,".csv"))
    for(i in celltypes){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL_loo[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL_loo[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
truevalues<-phenomeg[-c(14,20),celltypes]

estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs_loo[[j]]),-c(14,20)], coefCellType = IDOLcoefs_loo[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)

for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    cleanIDOL_loo[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    cleanIDOL_loo[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}



#Strict
allcordcleanstr<-allcordclean[,-c(which(allcordclean$CellType=="CD8T" & allcordclean$CellType2=="CD4T"),
                                  which(allcordclean$CellType=="NK" & allcordclean$CellType2=="CD8T"),
                                  which(allcordclean$CellType=="CD4T" & allcordclean$CellType2=="CD8T"))]
table(allcordcleanstr$Study, allcordcleanstr$CellType)
table(allcordcleanstr$CellType, allcordcleanstr$CellType2)
allcordcleanstrnoob<-preprocessNoob(allcordcleanstr)

Bakulski<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleanstrnoob[,allcordcleanstrnoob$Study=="Bakulski"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Bakulski$compTable[IDOLOptim517,3:9])
deGoede<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleanstrnoob[,allcordcleanstrnoob$Study=="deGoede"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(deGoede$compTable[IDOLOptim517,3:9])
Gervin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleanstrnoob[,allcordcleanstrnoob$Study=="Gervin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Gervin$compTable[IDOLOptim517,3:8])
Lin<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleanstrnoob[,allcordcleanstrnoob$Study=="Lin"], cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),probeSelect = "any") 
head(Lin$compTable[IDOLOptim517,3:8])
Allref<-FlowSorted.Blood.EPIC:::pickCompProbes(mSet = allcordcleanstrnoob, cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),probeSelect = "any") 
head(Allref$compTable[IDOLOptim517,3:9])

IDOLcoefs_strict<-list(Bakulski=as.matrix(Bakulski$compTable[IDOLOptim517,3:9]), deGoede=as.matrix(deGoede$compTable[IDOLOptim517,3:9]), Gervin=as.matrix(Gervin$compTable[IDOLOptim517,3:8]), 
                Lin=as.matrix(Lin$compTable[IDOLOptim517,3:8]), Allref=as.matrix(Allref$compTable[IDOLOptim517,3:9]))

cleanIDOL_strict<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                               "deGoede_R2", "deGoede_rmse",
                                                                               "Gervin_R2", "Gervin_rmse",
                                                                               "Lin_R2", "Lin_rmse",
                                                                               "Allref_R2", "Allref_rmse")))


for (j in Study[1:2]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs_strict[[j]]),-c(14,20)], coefCellType = IDOLcoefs_strict[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_strict_",j,".csv"))
    for(i in celltypes){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL_strict[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL_strict[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

for (j in Study[3:4]){
    
    
    truevalues<-phenomeg[-c(14,20),celltypes]
    
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs_strict[[j]]),-c(14,20)], coefCellType = IDOLcoefs_strict[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_strict_",j,".csv"))
    for(i in celltypes[1:6]){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL_strict[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL_strict[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
truevalues<-phenomeg[-c(14,20),celltypes]

estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(megdatanoob)[rownames(IDOLcoefs_strict[[j]]),-c(14,20)], coefCellType = IDOLcoefs_strict[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)
write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_strict_",j,".csv"))

for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    cleanIDOL_strict[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    cleanIDOL_strict[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}

save(IDOLcoefs, IDOLcoefs_loo, IDOLcoefs_strict, cleancoefs, uncleancoefs, Study, file = "Objects_for_IDOL_paper.RData")
save(cleanminfi, uncleanminfi, cleanIDOL, cleanIDOL_loo, cleanIDOL_strict, file="Results_R2_rmse.RData")




#QC ontrols
source("D:/OneDrive/OneDrive - Dartmouth College/Swimming pools/RICH/QCinfo2.R")

qcinfo<-QCinfo2(Referencecord)


load("D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive/IDOL cord/IDOL optimized DMR library_550.RData")
allcordclean<-allrefnew[, allrefnew$Cord=="Yes" & allrefnew$CellType!="MIX" & allrefnew$CellType!="Endothelial" & allrefnew$CellType!="Epithelial" & allrefnew$CellType!="Stromal"]
allcordcleannoob<-preprocessNoob(allcordclean)
IDOLcordcomp.Table<-IDOLcoefs[["Allref"]]
save(allcordclean, IDOLOptim517, allcordcleannoob, IDOLcordcomp.Table, file="D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive/IDOL cord/Reference_cord_blood_clean.RData")

allcordunclean<-allrefnew[, allrefnew$Cord=="Yes" & allrefnew$CellType!="Endothelial" & allrefnew$CellType!="Epithelial" & allrefnew$CellType!="Stromal"]
mixsamples<-colnames(allcordunclean[,allcordunclean$CellType=="MIX"])
mixedsamples<-as.data.frame(pData(allcordunclean[,mixsamples]))
swapsamples<-colnames(allcordunclean[,(allcordunclean$CellType=="CD8T" & allcordunclean$CellType2=="CD4T") |
                                         (allcordunclean$CellType=="CD4T" & allcordunclean$CellType2=="CD8T") | 
                                         (allcordunclean$CellType=="NK" & allcordunclean$CellType2=="CD8T") ])
swappedsamples<-as.data.frame(pData(allcordunclean[,swapsamples]))
save(allcordunclean, mixedsamples, swappedsamples, file="D:/OneDrive/OneDrive - Dartmouth College/Cord blood letter to the editor/archive/IDOL cord/Reference_cord_blood_unclean.RData")

