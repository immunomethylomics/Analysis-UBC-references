##################################
#
#Cord blood references deconvolution
#GenR validation of IDOL
#Author: Lucas A Salas
#Dartmouth College
#Date: 10/23/2018
###################################

#STEP 1:
#Load the RData containing the coefficient matrices

load("Objects_for_IDOL_paper.RData")

#STEP 2:
#Load the libraries required and your RgSet object here named "rgset" replace for yours
library(FlowSorted.Blood.EPIC)
library(minfi)
is(rgset, "RGChannelset")
#TRUE
noobdb<-preprocessNoob(rgset)
#Option if you have a Methylset object laready background corrected unmarke the line below and use it instead.
#noobdb<-yourobject
is(noobdb, "MethylSet")
#TRUE

#If FALSE please replace for the proper object.
#Do not use any additional normalization or filtering here

#STEP 3: select your true values (FACS information)
celltypes<-c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
truevalues<-phenotype[,celltypes]#replace phenotype for your object containing the FACS information, columns MUST have the same names as the celltypes for the comparisons
truevalues<-truevalues[order(rownames(truevalues)),]
write.csv(truevalues, file="True_values.csv")
#The rownames of your truevalues should be the same ID as your colnames id in your noobdb
noobdb<-noobdb[,order(colnames(noobdb))]
identical(rownames(truevalues), colnames(noobdb))
#TRUE
#If false correct the IDs before running the loop and corroborate the order of your samples
#IF THE ORDER OR IDS ARE INCORRECT THE LOOP WILL RUN!!!!!!!

#STEP 4 Run the loop
rmse<- function(error){
    sqrt(mean(error^2))
}

cleanminfi<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                "deGoede_R2", "deGoede_rmse",
                                                                                "Gervin_R2", "Gervin_rmse",
                                                                                "Lin_R2", "Lin_rmse",
                                                                                "Allref_R2", "Allref_rmse")))


for (j in Study[1:2]){
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(cleancoefs[[j]]),], coefCellType = cleancoefs[[j]])
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
   estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(cleancoefs[[j]]),], coefCellType = cleancoefs[[j]])
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
estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(cleancoefs[[j]]),], coefCellType = cleancoefs[[j]])
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


#Unclean
uncleanminfi<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                  "deGoede_R2", "deGoede_rmse",
                                                                                  "Gervin_R2", "Gervin_rmse",
                                                                                  "Lin_R2", "Lin_rmse",
                                                                                  "Allref_R2", "Allref_rmse")))

for (j in Study[1:2]){
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(uncleancoefs[[j]]),], coefCellType = uncleancoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("minfi_uncleaned_",j,".csv"))
    
    for(i in celltypes){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        uncleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        uncleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

for (j in Study[3:4]){
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(uncleancoefs[[j]]),], coefCellType = uncleancoefs[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("minfi_uncleaned_",j,".csv"))
    
    for(i in celltypes[1:6]){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        uncleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        uncleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(uncleancoefs[[j]]),], coefCellType = uncleancoefs[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)
write.csv(cellcount_EPICIDOL, file=paste0("minfi_uncleaned_",j,".csv"))

for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    uncleanminfi[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    uncleanminfi[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}


#IDOL

cleanIDOL<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                               "deGoede_R2", "deGoede_rmse",
                                                                               "Gervin_R2", "Gervin_rmse",
                                                                               "Lin_R2", "Lin_rmse",
                                                                               "Allref_R2", "Allref_rmse")))

for (j in Study[1:2]){
   estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs[[j]]),], coefCellType = IDOLcoefs[[j]])
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
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs[[j]]),], coefCellType = IDOLcoefs[[j]])
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
estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs[[j]]),], coefCellType = IDOLcoefs[[j]])
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


#IDOL leave one out (loo)
cleanIDOL_loo<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                   "deGoede_R2", "deGoede_rmse",
                                                                                   "Gervin_R2", "Gervin_rmse",
                                                                                   "Lin_R2", "Lin_rmse",
                                                                                   "Allref_R2", "Allref_rmse")))


for (j in Study[1:2]){
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs_loo[[j]]),], coefCellType = IDOLcoefs_loo[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_loo_",j,".csv"))
    for(i in celltypes){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL_loo[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL_loo[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

for (j in Study[3:4]){
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs_loo[[j]]),], coefCellType = IDOLcoefs_loo[[j]])
    cellcount_EPICIDOL<-estimatestmp*100
    rownames(cellcount_EPICIDOL)<-rownames(truevalues)
    write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_loo_",j,".csv"))
    for(i in celltypes[1:6]){
        lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
        rmse(lm.D9$residuals)
        summary(lm.D9)$r.squared
        cleanIDOL_loo[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
        cleanIDOL_loo[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
    }
}

j<-"Allref"
estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs_loo[[j]]),], coefCellType = IDOLcoefs_loo[[j]])
cellcount_EPICIDOL<-estimatestmp*100
rownames(cellcount_EPICIDOL)<-rownames(truevalues)
write.csv(cellcount_EPICIDOL, file=paste0("IDOL_cleaned_loo_",j,".csv"))
for(i in celltypes){
    lm.D9 <- lm(cellcount_EPICIDOL[,i]~truevalues[,i] )
    rmse(lm.D9$residuals)
    summary(lm.D9)$r.squared
    cleanIDOL_loo[i,paste0(j,"_R2")]=round(summary(lm.D9)$r.squared[1]*100,1)
    cleanIDOL_loo[i,paste0(j,"_rmse")]=round(rmse(lm.D9$residuals)[1],2)
}

#Strict
cleanIDOL_strict<-matrix(data = NA, ncol = 10, nrow = 7, dimnames = list(celltypes, c("Bakulski_R2", "Bakulski_rmse",
                                                                                      "deGoede_R2", "deGoede_rmse",
                                                                                      "Gervin_R2", "Gervin_rmse",
                                                                                      "Lin_R2", "Lin_rmse",
                                                                                      "Allref_R2", "Allref_rmse")))
for (j in Study[1:2]){
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs_strict[[j]]),], coefCellType = IDOLcoefs_strict[[j]])
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
    estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs_strict[[j]]),], coefCellType = IDOLcoefs_strict[[j]])
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
estimatestmp<-FlowSorted.Blood.EPIC:::projectCellType(Y = getBeta(noobdb)[rownames(IDOLcoefs_strict[[j]]),], coefCellType = IDOLcoefs_strict[[j]])
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

save(cleanminfi, uncleanminfi, cleanIDOL, cleanIDOL_loo, cleanIDOL_strict, file="Results_R2_rmse_GenR.RData")

#STEP 5 Please send the results all csv files and the "Results_R2_rmse_GenR.RData" to Kristina for further processing, flagging and QC filtering

Referencecord<-allcordclean
save(IDOLOptim517, Referencecord, file="Reference_cord_clean.RData")
