#XENOPUS PROJECT - MASTER CODE
#Last Updated: 8/15/19
#Aaron Ta

library(tidyr)
library("DESeq2")
library("tidyverse")
library("VennDiagram")

writeGeneList <- function(list, file) {
  write.table(list,file,quote=F,row.names = F,col.names = F,sep="\t")
}

###READ IN FEATURE COUNT FILES###
filenames <- list.files(path="/Users/Aaron/Desktop/STARAlignment4_10_19/10_17_19_9.2_exon_GTF/Counts/",
                        pattern="*.txt")
names <-substring(filenames, regexpr("_", filenames) + 1)
names <-sub("\\..*", "", names)
names<-c("MB44_I","MB44_II","MB44_III","FB46_I","FB46_II","FB46_III","HB46_I","HB46_II","HB46_III","MB46_I","MB46_II","MB46_III","SC46_I","SC46_II","SC46_III","FB49_I","FB49_II",
         "HB49_I","HB49_II","MB49_I","MB49_II","SC49_I","SC49_II","MB55_I","MB55_II","MB55_III","MB61_I","MB61_II","MB61_III","MB66_I")
Xennames<-names
for(i in 1:length(names)){
  assign(names[i], read.table(paste("/Users/Aaron/Desktop/STARAlignment4_10_19/10_17_19_9.2_exon_GTF/Counts/",filenames[i],sep=""),sep="\t"))
}


###CONVERSION TABLE CONSTRUCTION###

#Read in the X. laevis transcript To X. tropicalis transcript table (from Briggs et al.)
LaevTranscriptToTropTranscript<-read.csv("~/Desktop/XenConversionFiles/XeLaevAllBest.txt",sep="\t", header=FALSE)

#Read in the X. laevis GTF file and construct a X. laevis transcript to X. laevis symbol table
GTF<-read.table("~/Downloads/XENLA_9.2_Xenbase.GTF",sep="\t")
GTFsplit<-separate(GTF,V9,c("gene_id","gene_name","transcript_id","transcript_name"),sep="; ")[,10:11]
GTFsplit<-as.data.frame(unique(GTFsplit))
GTFsplit[,1]<-gsub("gene_name ","",GTFsplit[,1])
GTFsplit[,2]<-gsub("transcript_id ","",GTFsplit[,2])
LaevTranscriptToLaevSymbol<-cbind(GTFsplit[,2],GTFsplit[,1])

#Read in the X. tropicalis gene symbol to H. sapiens gene symbol table
TropToHumanGeneSymbol<-read.csv("~/Desktop/XenConversionFiles/XeLaevAllBestSplit/found/done/human/XeTropHumanGeneSymbol.txt", header=FALSE, sep="\t", fill=TRUE, blank.lines.skip = FALSE)
TropToHumanGeneSymbol<-as.matrix(TropToHumanGeneSymbol)
dim(TropToHumanGeneSymbol)
dim(LaevTranscriptToTropTranscript)

#Convert incorrect gene symbols in above table
TropToHumanGeneSymbol[TropToHumanGeneSymbol=="3-Sep"]<-"SEPT3"
TropToHumanGeneSymbol[TropToHumanGeneSymbol=="7-Sep"]<-"SEPT7"
TropToHumanGeneSymbol[TropToHumanGeneSymbol=="6-Mar"]<-"MARCH6"
TropToHumanGeneSymbol[TropToHumanGeneSymbol=="10-Mar"]<-"MARCH10"

#Make table connecting X. laevis transcript to X. tropicalis transcript to H. sapiens gene symbol
TropToHumanGeneSymbol<-cbind(LaevTranscriptToTropTranscript[,1],TropToHumanGeneSymbol[,1],LaevTranscriptToTropTranscript[,2:12],TropToHumanGeneSymbol[,2])

#Connect X. laevis symbol to X. laevis transcript (...to H. sapiens gene symbol by the previous table)
colnames(LaevTranscriptToLaevSymbol)<-cbind("LaevisGeneID","LaevisGeneSymbol")
colnames(TropToHumanGeneSymbol)<-cbind("LaevisGeneID2","TropGeneSymbol","TropGeneID","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","SapiensGeneSymbol")
LaevSymbolToSapiensGeneSymbol<-merge(x=LaevTranscriptToLaevSymbol,y=TropToHumanGeneSymbol,by.x="LaevisGeneID",by.y="LaevisGeneID2", all=T, stringsAsFactors = FALSE)
LaevSymbolToSapiensGeneSymbol$TropGeneSymbol<-as.character(LaevSymbolToSapiensGeneSymbol$TropGeneSymbol)
LaevSymbolToSapiensGeneSymbol$TropGeneID<-as.character(LaevSymbolToSapiensGeneSymbol$TropGeneID)
LaevSymbolToSapiensGeneSymbol$SapiensGeneSymbol<-as.character(LaevSymbolToSapiensGeneSymbol$SapiensGeneSymbol)

#Construct a table w/ no blanks
#Take care of NAs in dataset (retain all genes present in genome)
LaevToHumanGeneSymbolNoBlanks<-LaevSymbolToSapiensGeneSymbol
LaevToHumanGeneSymbolNoBlanks$TropGeneSymbol[is.na(LaevToHumanGeneSymbolNoBlanks$TropGeneSymbol)] <- ""
LaevToHumanGeneSymbolNoBlanks$TropGeneID[is.na(LaevToHumanGeneSymbolNoBlanks$TropGeneID)] <- ""
LaevToHumanGeneSymbolNoBlanks$SapiensGeneSymbol[is.na(LaevToHumanGeneSymbolNoBlanks$SapiensGeneSymbol)] <- ""
for (i in 5:14) {
  LaevToHumanGeneSymbolNoBlanks[is.na(LaevToHumanGeneSymbolNoBlanks[,i]),i] <- 0
}

#Fill in blanks with Trop Gene Symbol
for (row in 1:nrow(LaevToHumanGeneSymbolNoBlanks)) {
  if (LaevToHumanGeneSymbolNoBlanks[row,15]=="" & LaevToHumanGeneSymbolNoBlanks[row,3]!="") {
    LaevToHumanGeneSymbolNoBlanks[row,15]=as.character(LaevToHumanGeneSymbolNoBlanks[row,3])
  }
}

#Fill in blanks with Laevis Gene Symbol
for (row in 1:nrow(LaevToHumanGeneSymbolNoBlanks)) {
  if (LaevToHumanGeneSymbolNoBlanks[row,15]=="" & LaevToHumanGeneSymbolNoBlanks[row,2]!="") {
    LaevToHumanGeneSymbolNoBlanks[row,15]=as.character(LaevToHumanGeneSymbolNoBlanks[row,2])
  }
}

#Pull unique Laev Gene -> Trop Genes
uLaevToTrop<-LaevSymbolToSapiensGeneSymbol[order(LaevSymbolToSapiensGeneSymbol$bitscore,decreasing=TRUE),]
uLaevToTrop<-uLaevToTrop[order(uLaevToTrop$LaevisGeneSymbol, uLaevToTrop$evalue),]
uLaevToTrop<-uLaevToTrop[!duplicated(uLaevToTrop[,2]),]
LaevToTropRaw<-uLaevToTrop
LaevToTropRaw<-cbind(as.character(LaevToTropRaw$LaevisGeneSymbol),as.character(LaevToTropRaw$TropGeneSymbol))
colnames(LaevToTropRaw)<-cbind("LaevGeneSymbol","TropGeneSymbol")
LaevToTropRaw[which(is.na(LaevToTropRaw))]<-""
LaevToTrop<-LaevToTropRaw

#Pull unique Laev Gene->Human Genes
uLaevToHuman<-LaevToHumanGeneSymbolNoBlanks[order( LaevToHumanGeneSymbolNoBlanks$bitscore,decreasing=TRUE),]
uLaevToHuman<-uLaevToHuman[order(uLaevToHuman$LaevisGeneSymbol, uLaevToHuman$evalue),]
uLaevToHuman<-uLaevToHuman[!duplicated(uLaevToHuman[,2]),]
LaevToHumanGeneSymbolNoBlanks<-uLaevToHuman

#LaevTranscript To Human GeneSymbol, simplified w/ no blanks (filled in from other species)
LaevToHuman<-cbind(as.character(LaevToHumanGeneSymbolNoBlanks$LaevisGeneSymbol),as.character(LaevToHumanGeneSymbolNoBlanks$SapiensGeneSymbol))
colnames(LaevToHuman)<-cbind("LaevGeneSymbol","SapiensGeneSymbol")

#LaevTranscript To Human GeneSymbol, simplified w/ blanks (human gene symbols only)
uLaevToHumanRaw<-LaevSymbolToSapiensGeneSymbol[order( LaevSymbolToSapiensGeneSymbol$bitscore,decreasing=TRUE),]
uLaevToHumanRaw<-uLaevToHumanRaw[order(uLaevToHumanRaw$LaevisGeneSymbol, uLaevToHumanRaw$evalue),]
uLaevToHumanRaw<-uLaevToHumanRaw[!duplicated(uLaevToHumanRaw[,2]),]
LaevToHumanRaw<-uLaevToHumanRaw
LaevToHumanRaw<-cbind(as.character(LaevToHumanRaw$LaevisGeneSymbol),as.character(LaevToHumanRaw$SapiensGeneSymbol))
colnames(LaevToHumanRaw)<-cbind("LaevGeneSymbol","SapiensGeneSymbol")
LaevToHumanRaw[which(is.na(LaevToHumanRaw))]<-""

###DIFFERENTIAL EXPRESSION ANALYSIS - BY STAGE - IN MIDBRAIN###

wt_MBall<-vector()
#Extract columns w/ midbrain data
MBnames<-grep("MB",names)
for(j in MBnames){
  wt_MBall<-cbind(wt_MBall,eval(as.name(names[j]))[,2])
}
colnames(wt_MBall)<-names[MBnames]
rownames(wt_MBall) <- MB44_I[,1]
#Remove MB49_II for being an outlier
wt_MBall<-wt_MBall[,-which(names[MBnames]=="MB49_II")]

coldataMB <- c("stMB44", "stMB44", "stMB44", "stMB46", "stMB46", "stMB46", "stMB49", "stMB55", "stMB55", "stMB55", "stMB61", "stMB61", "stMB61", "stMB66")
coldataMB <- as.data.frame(coldataMB)
rownames(coldataMB)<-(names[MBnames])[-which(names[MBnames]=="MB49_II")]
colnames(coldataMB) <- "stage"

#DESeq for MB
ddsMB <- DESeqDataSetFromMatrix(countData = wt_MBall,
                              colData = coldataMB,
                              design = ~ stage)
keep <- rowSums(counts((ddsMB))) >= 14
ddsMB <- ddsMB[keep,]
ddsMB <- DESeq(ddsMB)

#DESeq for MB, no 49/66 for PCA
ddsMB_PCA <- DESeqDataSetFromMatrix(countData = wt_MBall[,c(-7,-14)],
                                colData = coldataMB[c(-7,-14),,drop=F],
                                design = ~ stage)
keep <- rowSums(counts((ddsMB_PCA))) >= 14
ddsMB_PCA <- ddsMB_PCA[keep,]
ddsMB_PCA <- DESeq(ddsMB_PCA)

#Two-stage comparison results
resMB44v46<-results(ddsMB,c("stage","stMB46","stMB44"))
resMB44v49<-results(ddsMB,c("stage","stMB49","stMB44"))
resMB44v55<-results(ddsMB,c("stage","stMB55","stMB44"))
resMB44v61<-results(ddsMB,c("stage","stMB61","stMB44"))
resMB44v66<-results(ddsMB,c("stage","stMB66","stMB44"))

resMB46v49<-results(ddsMB,c("stage","stMB49","stMB46"))
resMB46v55<-results(ddsMB,c("stage","stMB55","stMB46"))
resMB46v61<-results(ddsMB,c("stage","stMB61","stMB46"))
resMB46v66<-results(ddsMB,c("stage","stMB66","stMB46"))

resMB49v55<-results(ddsMB,c("stage","stMB55","stMB49"))
resMB49v61<-results(ddsMB,c("stage","stMB61","stMB49"))
resMB49v66<-results(ddsMB,c("stage","stMB66","stMB49"))

resMB55v61<-results(ddsMB,c("stage","stMB61","stMB55"))
resMB55v66<-results(ddsMB,c("stage","stMB66","stMB55"))

resMB61v66<-results(ddsMB,c("stage","stMB66","stMB61"))
resZD$log2FoldChange
StageComparisons<-c("resMB44v46","resMB44v49","resMB44v55","resMB44v61","resMB44v66",
                    "resMB46v49","resMB46v55","resMB46v61","resMB46v66",
                    "resMB49v55","resMB49v61","resMB61v66",
                    "resMB55v61","resMB61v66",
                    "resMB61v66")

#Store normalized counts for stage comparisons
wt_MBn<-counts(ddsMB, normalized=T)

#Set NA adjusted p-values to 1 (not significant)
for (k in 1:length(StageComparisons)) {
  tmpRes<-eval(as.name(StageComparisons[k]))
  tmpRes$padj[is.na(tmpRes$padj)] <- 1
  assign(StageComparisons[k],tmpRes) #Reassign modified table back into its original variable
}

#Make lists of significantly DE'd genes for each comparison: adjusted p-value < 0.05 and log2fold change magnitude > 2
#List stored in: resMB##v##genes (ie: resMB46v61genes)
sigStageGeneLists<-vector()
sigStageGeneListsUR<-vector() #significantly upregulated genes only
sigStageGeneListsDR<-vector() #significantly downregulated genes only
#By padj < 0.05, lfc > 2
for (l in 1:length(StageComparisons)) {
  tmpRes<-eval(as.name(StageComparisons[l]))
  sigStageGeneLists<-c(sigStageGeneLists,paste(StageComparisons[l],"genes",sep=""))
  sigStageGeneListsUR<-c(sigStageGeneListsUR,paste(StageComparisons[l],"genesUR",sep=""))
  sigStageGeneListsDR<-c(sigStageGeneListsDR,paste(StageComparisons[l],"genesDR",sep=""))
  assign(sigStageGeneLists[l], rownames(tmpRes[tmpRes$padj < 0.05 & abs(tmpRes$log2FoldChange) > 2,]))
  assign(sigStageGeneListsUR[l], rownames(tmpRes[tmpRes$padj < 0.05 & tmpRes$log2FoldChange > 2,]))
  assign(sigStageGeneListsDR[l], rownames(tmpRes[tmpRes$padj < 0.05 & tmpRes$log2FoldChange < -2,]))
}

#By p > 0.0001 (currently using padj<0.05 + | l2fc |> 2, not this cutoff)
#for (l in 1:length(StageComparisons)) {
#  tmpRes<-eval(as.name(StageComparisons[l]))
#  sigStageGeneLists<-c(sigStageGeneLists,paste(StageComparisons[l],"genes",sep=""))
#  sigStageGeneListsUR<-c(sigStageGeneListsUR,paste(StageComparisons[l],"genesUR",sep=""))
#  sigStageGeneListsDR<-c(sigStageGeneListsDR,paste(StageComparisons[l],"genesDR",sep=""))
#  assign(sigStageGeneLists[l], rownames(tmpRes[tmpRes$padj < 0.0001 & abs(tmpRes$log2FoldChange) > 2,]),]))
#  assign(sigStageGeneListsUR[l], rownames(tmpRes[tmpRes$padj < 0.0001 & tmpRes$log2FoldChange > 0,]))
#  assign(sigStageGeneListsDR[l], rownames(tmpRes[tmpRes$padj < 0.0001 & tmpRes$log2FoldChange > 0,]))
#}

#H. sapiens conversion for GO
#List stored in: resMB##v##genesHum (ie: resMB46v61genesHum)
sigStageGeneListsHum<-vector()
for (l in 1:length(sigStageGeneLists)) {
  tmpRes<-eval(as.name(sigStageGeneLists[l]))
  tmpRes<-unique(merge(tmpRes,LaevToHumanRaw,by.x=1,by.y="LaevGeneSymbol")[,2])
  sigStageGeneListsHum<-c(sigStageGeneListsHum,paste(sigStageGeneLists[l],"Hum",sep=""))
  assign(sigStageGeneListsHum[l], tmpRes)
}


#List of all genes that are significantly DE'd between (st44/st46) and (st55/st61)
resBGMBgenes<-rownames(resMB44v55)
resAllMBgenes<-union(resMB44v55genes,union(resMB44v61genes,union(resMB46v55genes,union(resMB44v46genes,union(resMB55v61genes,resMB46v61genes)))))
resAllMBgenesHum<-union(resMB44v55genesHum,union(resMB44v61genesHum,union(resMB46v55genesHum,union(resMB44v46genesHum,union(resMB55v61genesHum,resMB46v61genesHum)))))

###DIFFERENTIAL EXPRESSION ANALYSIS - BY REGION - IN STAGE 46###

wt_46all<-vector()
#Extract columns w/ st46 data
st46names<-grep("46",names)
for(j in st46names){
  wt_46all<-cbind(wt_46all,eval(as.name(names[j]))[,2])
}
colnames(wt_46all)<-names[st46names]
rownames(wt_46all) <- MB46_I[,1]

coldata46 <- c("stFB46", "stFB46", "stFB46", "stHB46", "stHB46", "stHB46", "stMB46", "stMB46", "stMB46","stSC46","stSC46","stSC46")
coldata46 <- as.data.frame(coldata46)
rownames(coldata46)<-(names[st46names])
colnames(coldata46) <- "region"

#DESeq for region
dds46 <- DESeqDataSetFromMatrix(countData = wt_46all,
                                colData = coldata46,
                                design = ~ region)
keep <- rowSums(counts((dds46))) >= 10
dds46 <- dds46[keep,]
dds46 <- DESeq(dds46)

#Two-region comparison results
res46FBvMB<-results(dds46,c("region","stFB46","stMB46"))
res46FBvHB<-results(dds46,c("region","stFB46","stHB46"))
res46MBvHB<-results(dds46,c("region","stMB46","stHB46"))
res46FBvSC<-results(dds46,c("region","stFB46","stSC46"))
res46MBvSC<-results(dds46,c("region","stMB46","stSC46"))
res46HBvSC<-results(dds46,c("region","stHB46","stSC46"))

RegionComparisons<-c("res46FBvMB","res46FBvHB","res46MBvHB","res46FBvSC","res46MBvSC","res46HBvSC")

#Store normalized counts for region comparison
wt_46n<-counts(dds46, normalized=T)

#Set NA adjusted p-values to 1 (not significant)
for (k in 1:length(RegionComparisons)) {
  tmpRes<-eval(as.name(RegionComparisons[k]))
  tmpRes$padj[is.na(tmpRes$padj)] <- 1
  assign(RegionComparisons[k],tmpRes) #Reassign modified table back into its original variable
}

#Make lists of significantly DE'd genes for each comparison: adjusted p-value < 0.05 and log2fold change magnitude > 2
#List stored in: res46##v##genes (ie: res46FBvMBgenes)
sigRegionGeneLists<-vector()
sigRegionGeneListsUR<-vector() #significantly upregulated genes only
sigRegionGeneListsDR<-vector() #significantly downregulated genes only
for (l in 1:length(RegionComparisons)) {
  tmpRes<-eval(as.name(RegionComparisons[l]))
  sigRegionGeneLists<-c(sigRegionGeneLists,paste(RegionComparisons[l],"genes",sep=""))
  sigRegionGeneListsUR<-c(sigRegionGeneListsUR,paste(RegionComparisons[l],"genesUR",sep=""))
  sigRegionGeneListsDR<-c(sigRegionGeneListsDR,paste(RegionComparisons[l],"genesDR",sep=""))
  assign(sigRegionGeneLists[l], rownames(tmpRes[tmpRes$padj < 0.05 & abs(tmpRes$log2FoldChange) > 2,]))
  assign(sigRegionGeneListsUR[l], rownames(tmpRes[tmpRes$padj < 0.05 & tmpRes$log2FoldChange > 2,]))
  assign(sigRegionGeneListsDR[l], rownames(tmpRes[tmpRes$padj < 0.05 & tmpRes$log2FoldChange < -2,]))
}

#List of genes that are significantly DE'd in at least one comparison
resBG46genes<-rownames(res46FBvHB)
resAll46genes<-union(res46FBvMBgenes,union(res46FBvHBgenes,union(res46MBvHBgenes,union(res46FBvSCgenes,union(res46MBvSCgenes,res46HBvSCgenes)))))

res46FBgenes<-intersect(res46FBvMBgenes,intersect(res46FBvHBgenes,res46FBvSCgenes))
res46MBgenes<-intersect(res46FBvMBgenes,intersect(res46MBvHBgenes,res46MBvSCgenes))
res46HBgenes<-intersect(res46FBvHBgenes,intersect(res46MBvHBgenes,res46HBvSCgenes))
res46SCgenes<-intersect(res46FBvSCgenes,intersect(res46MBvSCgenes,res46HBvSCgenes))

###CLUSTER PROFILING - BY STAGE IN MB###
designMB <- as.data.frame(colData(ddsMB))
clustersMB <- DEGreport::degPatterns(wt_MBn[rownames(wt_MBn) %in% resAllMBgenes,], metadata = designMB, time = "stage", col=NULL)
clusterMB_groups <- clustersMB$df
#Clustering after taking the log2 of the DESeq2-normalized counts
#clustersMBlog2 <- DEGreport::degPatterns(log2(wt_MBn[rownames(wt_MBn) %in% resAllMBgenes,]), metadata = designMB, time = "stage", col=NULL)

###CLUSTER PROFILING - BY STAGE IN MB (NO 49, 66)###
designMB <- as.data.frame(colData(ddsMB)[c(-7,-14),])
clustersMB <- DEGreport::degPatterns(wt_MBn[rownames(wt_MBn) %in% resAllMBgenes,c(-7,-14)], metadata = designMB, time = "stage", col=NULL, reduce = TRUE)
clusterMB_groups <- clustersMB$df

#Sample of how to extract a group from above (in this case, cluster 1)
groupMB1 <- clusterMB_groups %>%
  filter(cluster == 1)

#Make gene lists for each group in the format "groupMB#" (ie: groupMB1)
groupMBLists<-vector()
for(i in unique(clusterMB_groups$cluster)){
  groupMBLists<-c(groupMBLists,(paste("groupMB",i,sep="")))
}

for (l in 1:length(groupMBLists)) {
  assign(groupMBLists[l], clusterMB_groups %>% filter(cluster == l))
  assign(groupMBLists[l], eval(as.name(groupMBLists[l]))[,1])
}

wt_All<-vector()
#Extract columns w/ midbrain data
Allnames<-names
for(jj in Allnames){
  wt_All<-cbind(wt_All,eval(as.name(jj))[,2])
}
colnames(wt_All)<-Allnames
rownames(wt_All) <- MB44_I[,1]

coldataAll <- c("stMB44", "stMB44", "stMB44", "stFB46", "stFB46", "stFB46", "stHB46", "stHB46", "stHB46", "stMB46", "stMB46", "stMB46", "stSC46", "stSC46", "stSC46", "stFB49", "stFB49", "stHB49", "stHB49", "stMB49", "stMB49", "stSC49", "stSC49", "stMB55", "stMB55", "stMB55", "stMB61", "stMB61", "stMB61", "stMB66")
coldataAll <- as.data.frame(coldataAll)
colnames(coldataAll) <- "condition"

#DESeq for All
ddsAll <- DESeqDataSetFromMatrix(countData = wt_All,
                                colData = coldataAll,
                                design = ~ condition)
keep <- rowSums(counts((ddsAll))) >= 10
ddsAll <- ddsAll[keep,]
ddsAll <- DESeq(ddsAll)

vsd <- vst(ddsAll, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"),ntop = 1000)

vsd <- vst(ddsMB_PCA, blind=FALSE)
plotPCA(vsd, intgroup=c("stage"),ntop = 1000)

vsd <- vst(dds46, blind=FALSE)
plotPCA(vsd, intgroup=c("region"),ntop = 1000)

###CLUSTER PROFILING - BY REGION IN st46###
design46 <- as.data.frame(colData(dds46))
clusters46 <- DEGreport::degPatterns(wt_46n[rownames(wt_46n) %in% resAll46genes,], metadata = design46, time = "region", col=NULL)
cluster46_groups<-clusters46$df
#Clustering after taking the log2 of the DESeq2-normalized counts
#cluster46log2 <- DEGreport::degPatterns(log2(wt_46n[rownames(wt_46n) %in% resAll46genes,]), metadata = design46, time = "stage", col=NULL)

#Sample of how to extract a group from above (in this case, cluster 1)
group46_1 <- cluster46_groups %>%
  filter(cluster == 1)

#Make gene lists for each group in the format "group46_#" (ie: group46_1)
group46Lists<-vector()
for(i in unique(cluster46_groups$cluster)){
  group46Lists<-c(group46Lists,(paste("group46_",i,sep="")))
}

for (l in 1:length(group46Lists)) {
  assign(group46Lists[l], cluster46_groups %>% filter(cluster == l))
  assign(group46Lists[l], eval(as.name(group46Lists[l]))[,1])
}

#Make mean, SD tables for graphs
wt_MBn<-counts(ddsMB, normalized=T)
wt_MBMean<-cbind(rowMeans(wt_MBn[,1:3]),rowMeans(wt_MBn[,4:6]),wt_MBn[,7],rowMeans(wt_MBn[,8:10]),rowMeans(wt_MBn[,11:13]),wt_MBn[,14])
colnames(wt_MBMean)<-c("MB44","MB46","MB49","MB55","MB61","MB66")
wt_MBSD<-cbind(rowSds(wt_MBn[,1:3]),rowSds(wt_MBn[,4:6]),0,rowSds(wt_MBn[,8:10]),rowSds(wt_MBn[,11:13]),0)
rownames(wt_MBSD)<-rownames(wt_MBMean)
colnames(wt_MBSD)<-c("MB44","MB46","MB49","MB55","MB61","MB66")


wt_46Mean<-cbind(rowMeans(wt_46n[,1:3]),rowMeans(wt_46n[,4:6]),rowMeans(wt_46n[,7:9]),rowMeans(wt_46n[,10:12]))
colnames(wt_46Mean)<-c("FB","HB","MB","SC")
wt_46SD<-cbind(rowMeans(wt_46n[,1:3]),rowMeans(wt_46n[,4:6]),rowMeans(wt_46n[,7:9]),rowMeans(wt_46n[,10:12]))
colnames(wt_46SD)<-c("FB","HB","MB","SC")

###DIFFERENTIAL EXPRESSION ANALYSIS - BY STAGE - st46v49###

wt_4649<-vector()
#Extract columns w/ midbrain data
names4649<-as.integer(c(grep("46",names),grep("49",names)))
for(j in names4649){
  wt_4649<-cbind(wt_4649,eval(as.name(names[j]))[,2])
}
colnames(wt_4649)<-names[names4649]
rownames(wt_4649) <- MB44_I[,1]
#Remove MB49_II for being an outlier
#wt_MBall<-wt_MBall[,-which(names[MBnames]=="MB49_II")]

coldata4649 <- c("46", "46", "46", "46", "46", "46", "46", "46", "46", "46", "46", "46", "49", "49", "49", "49", "49", "49", "49", "49")
coldata4649 <- as.data.frame(coldata4649)
#rownames(coldata4649)<-(names[names4649])[-which(names[names4649]=="MB49_II")]
colnames(coldata4649) <- "stage"

#DESeq for MB
dds4649 <- DESeqDataSetFromMatrix(countData = wt_4649,
                                colData = coldata4649,
                                design = ~ stage)
keep <- rowSums(counts((dds4649))) >= 10
dds4649 <- dds4649[keep,]
dds4649 <- DESeq(dds4649)

res46v49<-results(dds4649,c("stage","46","49"))

#Store normalized counts for region comparison
wt_46v49n<-counts(dds4649, normalized=T)

#Set NA adjusted p-values to 1 (not significant)
res46v49$padj[is.na(res46v49$padj)] <- 1

#Make list of significantly DE'd genes for each comparison: adjusted p-value < 0.05 and log2fold change magnitude > 2
res46v49genes<-rownames(res46v49[res46v49$padj < 0.05 & abs(res46v49$log2FoldChange) > 2,])

######

rv <- rowVars(assay(ddsMB))
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(dds)[1:1000,]))


#Save plots as tiff
  #tiff(name, height = h, width = w, units='cm', compression = "lzw", res = 300)
#Save plots as eps
  #setEPS()
  #postscript(name, height = h, width = w)

