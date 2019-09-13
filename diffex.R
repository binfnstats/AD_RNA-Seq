library("VennDiagram")
library("RUVSeq")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("topGO")
library("sleuth")
library("biomaRt")
library(WGCNA)
library("ggbiplot")
library(reshape)
library(gplots)
library(ops)
library(calibrate)
library(biomaRt)
library(sva)
library(ggplot2)
library("corrplot")
library(gage)
library("pathview")

covar=read.delim("../pd_rna_covar.txt")
View(covar)
covar=covar[order(covar$region,covar$condition),]

reads=read.delim("pd_rna_rsem_gene.tsv",row.names = 1)
reads_1=reads[,which(!colSums(reads)<200000)]

reads_wo_ctrl=reads[,-1]
reads_wo_ctrl=reads_wo_ctrl[,covar$ALFRED]
colnames(reads_wo_ctrl)=covar$sample
reads_wo_ctrl1=reads_wo_ctrl[,match(covar$sample,colnames(reads_wo_ctrl))]
head(reads_wo_ctrl)
keep_soc_pd=apply(reads_wo_ctrl[,which(grepl("SOC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=1])>=5)
keep_soc_ctrl=apply(reads_wo_ctrl[,which(grepl("CTRL_SOC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=1])>=5)
keep_sfc_ctrl=apply(reads_wo_ctrl[,which(grepl("CTRL_SFC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=1])>=5)
keep_sfc_pd=apply(reads_wo_ctrl[,which(grepl("PD_SFC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=1])>=5)
keep_bf_pd=apply(reads_wo_ctrl[,which(grepl("PD_BF",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=1])>=5)
keep_bf_ctrl=apply(reads_wo_ctrl[,which(grepl("CTRL_BF",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=1])>=5)
keep=keep_bf_ctrl|keep_bf_pd|keep_sfc_ctrl|keep_sfc_pd|keep_soc_ctrl|keep_soc_pd
reads_fil=reads_wo_ctrl[keep,]
dim(reads_wo_ctrl)
dim(reads_fil)


good=goodGenes(datExpr0)
reads_good=datExpr0[,good]

y=DGEList(counts=t(reads_good)
	
y=DGEList(counts=reads_fil)
y=calcNormFactors(y,method = "TMM")

reads_norm<-cpm(y, normalized.lib.sizes=T)
datExpr0=t(reads_norm)
write.table(reads_norm,file="Normalised_reads.tsv",sep="\t")

####unsupervised clustering to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
dev.off()

###PCA analaysis
rsem.pca <- prcomp(t(reads_norm), scale = TRUE)

pdf(file = "PCA_rsem_rin.pdf")
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1,groups = covar$rin_cat, ellipse = TRUE, circle = TRUE, var.axes = T,labels=covar$No)
scale_color_discrete(name = '') +
theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

####Batch removal
design <- model.matrix(~1,data=covar$Condition)
dat1=as.matrix(reads_fil) 
###or 
dat1=reads_norm
combat <- ComBat(dat=dat1,batch=covar$Batch,mod=design,par.prior=TRUE)

####ANOVA analysis
z=reads_fil
m=melt(z)
colnames(m) <- c("sample_ID","counts")

dis <- rep(as.numeric(covar$Disease_merged, each=nrow(z)))
matrix <- data.frame(m,Disease=dis,ORF=orf,Region=reg,Patient=pat,Age=age,Gender=sex,PMD=pmi,RIN=rin)
fit1 <- lm(counts ~ Disease + ORF + Region + Patient + Age + Gender + PMD + RIN , data=matrix)
a <- anova(fit1)
nfac <- length(a[,1])-1
maxval = 100
pdf(file="Anovar.pdf")
nfac <- length(a[,1])-1
barplot(100*a$"Sum Sq"[1:nfac]/sum(a$"Sum Sq"[1:nfac]),names.arg=rownames(a[1:nfac,]),ylim=c(0,maxval),las=3)
dev.off()

design=model.matrix(~0+covar$partner,data=covar$condition)
y=estimateGLMCommonDisp(y,design,verbose=T)
y=estimateGLMTagwiseDisp(y,design)
fit=glmFit(y,design )
pdf(file="BCV.pdf")
plotBCV(y)
dev.off()
ctrl_vs_pd_lrt=glmLRT(fit,contrast=c(-1,-1,-1,1,1,1))
design
ctrl_bf_vs_pd_bf_lrt=glmLRT(fit,contrast=c(-1,0,0,1,0,0))
ctrl_sfc_vs_pd_sfc_lrt=glmLRT(fit,contrast=c(0,-1,0,0,1,0))
ctrl_soc_vs_pd_soc_lrt=glmLRT(fit,contrast=c(0,0,-1,0,0,1))


mart=useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "mar2016.archive.ensembl.org")

kg.hsa=kegg.gsets("hsa")
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.gs.sym <- lapply(kegg.gs, eg2sym)
colors=brewer.pal(8,"Set3")
x=covar$partner


x= c(rep("HC", length(which(grepl("HC",colnames(reads_fil))))), rep("ILB", length(which(grepl("ILB",colnames(reads_fil))))), rep("PD", length(which(grepl("PD",colnames(reads_fil))))))


keep_soc_pd=apply(reads_wo_ctrl[,which(grepl("PD_SOC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep_soc_ctrl=apply(reads_wo_ctrl[,which(grepl("CTRL_SOC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep_sfc_ctrl=apply(reads_wo_ctrl[,which(grepl("CTRL_SFC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep_sfc_pd=apply(reads_wo_ctrl[,which(grepl("PD_SFC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep_bf_pd=apply(reads_wo_ctrl[,which(grepl("PD_BF",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep_bf_ctrl=apply(reads_wo_ctrl[,which(grepl("CTRL_BF",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep=keep_bf_ctrl|keep_bf_pd|keep_sfc_ctrl|keep_sfc_pd|keep_soc_ctrl|keep_soc_pd
reads_fil=reads_wo_ctrl[keep,]
dim(reads_wo_ctrl)
dim(reads_fil)



pd_rna_rsem_gene_cpm=cpm(y, normalized.lib.sizes=T)
pd_rna_rsem_gene_cpm_t=t(cpm(y, normalized.lib.sizes=T))
dput(pd_rna_rsem_gene_cpm_t,file="harvard_kalli_iso_cpm_t.dput")
dput(pd_rna_rsem_gene_cpm,file="harvard_kalli_iso_cpm.dput")
z=x$label2
dput(z,file="pd_rna_label.dput")


ctrl_vs_pd_lrt=glmLRT(fit,contrast=c(1,1,1,-1,-1,-1))
ctrl_bf_vs_pd_bf_lrt=glmLRT(fit,contrast=c(1,0,0,-1,0,0))
ctrl_sfc_vs_pd_sfc_lrt=glmLRT(fit,contrast=c(0,1,0,0,-1,0))
ctrl_soc_vs_pd_soc_lrt=glmLRT(fit,contrast=c(0,0,1,0,0,-1))





keep_hc=apply(reads_wo_ctrl[,which(grepl("HC",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=20)
keep_ilb=apply(reads_wo_ctrl[,which(grepl("ILB",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=5)
keep_pd=apply(reads_wo_ctrl[,which(grepl("PD",colnames(reads_wo_ctrl)))],1,function(x) length(x[x>=10])>=2)
keep=keep_hc|keep_ilb|keep_pd
reads_fil=reads_wo_ctrl[keep,]
dim(reads_wo_ctrl)
dim(reads_fil)


for i in `cat test1`; do cp ./cki/p01/${i}_DEG.tsv results/p01/cki_${i}_DEG.tsv; done

results_DEG_mart <- getBM(attributes = c("ensembl_gene_id", "refseq_mrna","hgnc_symbol", "entrezgene", "description", "chromosome_name", "start_position", "end_position", "strand","go_id"), filters = "ensembl_gene_id", values = sub("\\..+","",rownames(DEG_pd_bf_vs_ctrl_bf_lrt)),mart = mart)

annotation= c(rep("HC", length(which(grepl("HC",rownames(datExpr))))), rep("ILB", length(which(grepl("ILB",rownames(datExpr))))), rep("PD", length(which(grepl("PD",rownames(datExpr))))))

annotation<-data.frame(covar$region,covar$condition,covar$sample)
colnames(annotation)<-c("Region","Condition","Sample")
rownames(annotation)<-covar$com
annotation<-subset(annotation,select=c(Region,Condition))

Condition = c("#e41a1c","#377eb8","#4daf4a")
names(Condition) = c("T","Sporadic","Orf")
Region = c("#f0f0f0", "#000000")
names(Region) = c("acute","chronic")

ann_col<-list("Condition"=Condition,"Region"=Region)

###Genes
for i in `cat test` ; do echo "

summary(dt_${i} <- decideTestsDGE(${i}, p=0.05, adjust=\"BH\"))
isDE_${i} <- as.logical(dt_${i})
DE_${i} <- rownames(y)[isDE_${i}]
${i}\$table\$FDR <- p.adjust(${i}\$table\$PValue, method=\"BH\")
DEG_${i} <- ${i}\$table[${i}\$table\$FDR < 0.05,]
DEG_${i}=DEG_${i}[DEG_${i}\$FDR > 0,]
DEG_${i}=DEG_${i}[DEG_${i}\$PValue < 0.05,]
#DEG_${i}=DEG_${i}[(DEG_${i}\$logFC<(-0.5))|(DEG_${i}\$logFC>0.5),]

results_DEG_${i} <- getBM(attributes = c(\"ensembl_gene_id\", \"hgnc_symbol\", \"entrezgene\", \"gene_biotype\",\"transcript_biotype\",\"description\", \"chromosome_name\", \"start_position\", \"end_position\", \"strand\",\"go_id\"), filters = \"ensembl_gene_id\", values = sub(\"\\\\..+\",\"\",rownames(DEG_${i})),mart = mart)

idx_DEG_${i}=match(sub(\"\\\\..+\",\"\",rownames(DEG_${i})),results_DEG_${i}\$ensembl_gene_id)

DEG_${i}\$ensembl_gene_id=results_DEG_${i}\$ensembl_gene_id[idx_DEG_${i}]
DEG_${i}\$ensembl_transcript_id=results_DEG_${i}\$ensembl_transcript_id[idx_DEG_${i}]
DEG_${i}\$hgnc_symbol=results_DEG_${i}\$hgnc_symbol[idx_DEG_${i}]
DEG_${i}\$entrezgene=results_DEG_${i}\$entrezgene[idx_DEG_${i}]
DEG_${i}\$gene_biotype=results_DEG_${i}\$gene_biotype[idx_DEG_${i}]
DEG_${i}\$transcript_biotype=results_DEG_${i}\$transcript_biotype[idx_DEG_${i}]
DEG_${i}\$description=results_DEG_${i}\$description[idx_DEG_${i}]
DEG_${i}\$chromosome_name=results_DEG_${i}\$chromosome_name[idx_DEG_${i}]
DEG_${i}\$start_position=results_DEG_${i}\$start_position[idx_DEG_${i}]
DEG_${i}\$end_position=results_DEG_${i}\$end_position[idx_DEG_${i}]
DEG_${i}\$strand=results_DEG_${i}\$strand[idx_DEG_${i}]
####DEG_${i}\$go_id=results_DEG_${i}\$go_id[idx_DEG_${i}]


write.table(DEG_${i}, file=\"${i}_DEG.tsv\", sep=\"\\t\")
top_${i}=rownames(DEG_${i})
${i}_cpm <- cpm(y)[top_${i}, ]
z_score_${i} <- ((${i}_cpm - rowMeans(${i}_cpm))/apply(${i}_cpm,1,sd))

write.table(${i}_cpm, file=\"${i}_DEG_cpm.tsv\", sep=\"\\t\", row.names =T)

unique_DEG_${i}=unique(DEG_${i}\$entrezgene)


pdf(\"${i}_heatmap.pdf\")
pheatmap(z_score_${i}, scale=\"row\", cluster_cols=T, cluster_rows=T, show_rownames=F,color = colorRampPalette(rev(brewer.pal(n = 11, name = \"RdYlBu\")))(100))
dev.off()

pdf(\"${i}_MA_plot.pdf\")
plotSmear(${i}, de.tags=DE_${i})
abline(h=c(-1,1), col=\"blue\")
dev.off()

pdf(\"${i}_pvalues.pdf\")
hist(${i}\$table\$PValue, breaks=seq(0, 1, 0.05))
dev.off()

pdf(\"${i}_lrt_RLE.pdf\")
plotRLE(${i}\$fitted.values, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

anno_bp <- annFUN.org(\"BP\", mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
allGenes_bp <- unique(unlist(anno_bp))
uniqueDEG_${i} <- unique(results_DEG_${i}\$ensembl_gene_id)
geneList_${i} <- factor(as.integer(allGenes_bp %in% uniqueDEG_${i}))
names(geneList_${i}) <- allGenes_bp

GOdata_bp_${i} <- new(\"topGOdata\", ontology = \"BP\", allGenes = geneList_${i}, nodeSize = 10, annot = annFUN.org, mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
resultFisher_bp_${i} <- runTest(GOdata_bp_${i}, algorithm = \"classic\", statistic = \"fisher\")
resultKS_bp_${i} <- runTest(GOdata_bp_${i}, algorithm = \"classic\", statistic = \"ks\")
resultKS.elim_bp_${i} <- runTest(GOdata_bp_${i}, algorithm = \"elim\", statistic = \"ks\")

pdf(\"${i}_gotree_BP.pdf\")
showSigOfNodes(GOdata_bp_${i}, score(resultFisher_bp_${i}), firstSigNodes = 10, useInfo = \"all\")
dev.off()

allRes_BP_${i} <- GenTable(GOdata_bp_${i}, classicFisher = resultFisher_bp_${i}, classicKS = resultKS_bp_${i}, elimKS = resultKS.elim_bp_${i}, orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 30)
for (i in 1:nrow(allRes_BP_${i})){ allRes_BP_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(allRes_BP_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(allRes_BP_${i})){ allRes_BP_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(allRes_BP_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(allRes_BP_${i}, file=\"${i}_toptree_BP.tsv\", sep=\"\\t\", row.names =F)

anno_cc <- annFUN.org(\"CC\", mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
allGenes_cc <- unique(unlist(anno_cc))
uniqueDEG_${i} <- unique(results_DEG_${i}\$ensembl_gene_id)
geneList_${i} <- factor(as.integer(allGenes_cc %in% uniqueDEG_${i}))
names(geneList_${i}) <- allGenes_cc

GOdata_cc_${i} <- new(\"topGOdata\", ontology = \"CC\", allGenes = geneList_${i}, nodeSize = 10, annot = annFUN.org, mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
resultFisher_cc_${i} <- runTest(GOdata_cc_${i}, algorithm = \"classic\", statistic = \"fisher\")
resultKS_cc_${i} <- runTest(GOdata_cc_${i}, algorithm = \"classic\", statistic = \"ks\")
resultKS.elim_cc_${i} <- runTest(GOdata_cc_${i}, algorithm = \"elim\", statistic = \"ks\")

pdf(\"${i}_gotree_CC.pdf\")
showSigOfNodes(GOdata_cc_${i}, score(resultFisher_cc_${i}), firstSigNodes = 10, useInfo = \"all\")
dev.off()

allRes_cc_${i} <- GenTable(GOdata_cc_${i}, classicFisher = resultFisher_cc_${i}, classicKS = resultKS_cc_${i}, elimKS = resultKS.elim_cc_${i}, orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 30)
for (i in 1:nrow(allRes_cc_${i})){ allRes_cc_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(allRes_cc_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(allRes_cc_${i})){ allRes_cc_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(allRes_cc_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(allRes_cc_${i}, file=\"${i}_toptree_CC.tsv\", sep=\"\\t\", row.names =F)

anno_mf <- annFUN.org(\"MF\", mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
allGenes_mf <- unique(unlist(anno_mf))
uniqueDEG_${i} <- unique(results_DEG_${i}\$ensembl_gene_id)
geneList_${i} <- factor(as.integer(allGenes_mf %in% uniqueDEG_${i}))
names(geneList_${i}) <- allGenes_mf

GOdata_mf_${i} <- new(\"topGOdata\", ontology = \"MF\", allGenes = geneList_${i}, nodeSize = 10, annot = annFUN.org, mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
resultFisher_mf_${i} <- runTest(GOdata_mf_${i}, algorithm = \"classic\", statistic = \"fisher\")
resultKS_mf_${i} <- runTest(GOdata_mf_${i}, algorithm = \"classic\", statistic = \"ks\")
resultKS.elim_mf_${i} <- runTest(GOdata_mf_${i}, algorithm = \"elim\", statistic = \"ks\")

pdf(\"${i}_gotree_MF.pdf\")
showSigOfNodes(GOdata_mf_${i}, score(resultFisher_mf_${i}), firstSigNodes = 10, useInfo = \"all\")
dev.off()

allRes_mf_${i} <- GenTable(GOdata_mf_${i}, classicFisher = resultFisher_mf_${i}, classicKS = resultKS_mf_${i}, elimKS = resultKS.elim_mf_${i}, orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 30)
for (i in 1:nrow(allRes_mf_${i})){ allRes_mf_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(allRes_mf_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(allRes_mf_${i})){ allRes_mf_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(allRes_mf_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(allRes_mf_${i}, file=\"${i}_toptree_MF.tsv\", sep=\"\\t\", row.names =F)

unique_DEG_${i}=unique(DEG_${i}\$entrezgene)
go_DEG_${i}=goana(unique_DEG_${i})

top_BP_DEG_${i}=topGO(go_DEG_${i},ont=\"BP\",n=30)
for (i in 1:nrow(top_BP_DEG_${i})){ top_BP_DEG_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(row.names(top_BP_DEG_${i})[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(top_BP_DEG_${i})){ top_BP_DEG_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(row.names(top_BP_DEG_${i})[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(top_BP_DEG_${i}, file=\"${i}_top1_BP.tsv\", sep=\"\\t\")

top_MF_DEG_${i}=topGO(go_DEG_${i},ont=\"MF\",n=30)
for (i in 1:nrow(top_MF_DEG_${i})){ top_MF_DEG_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(row.names(top_MF_DEG_${i})[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(top_MF_DEG_${i})){ top_MF_DEG_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(row.names(top_MF_DEG_${i})[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(top_MF_DEG_${i}, file=\"${i}_top1_MF.tsv\", sep=\"\\t\")

top_CC_DEG_${i}=topGO(go_DEG_${i},ont=\"CC\",n=30)
for (i in 1:nrow(top_CC_DEG_${i})){ top_CC_DEG_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(row.names(top_CC_DEG_${i})[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(top_CC_DEG_${i})){ top_CC_DEG_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(row.names(top_CC_DEG_${i})[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(top_CC_DEG_${i}, file=\"${i}_top1_CC.tsv\", sep=\"\\t\")



all_ensg_mart <- getBM(attributes = c(\"ensembl_gene_id\", \"entrezgene\"), filters = \"ensembl_gene_id\", values = sub(\"\\\\..+\",\"\",rownames(${i})),mart = mart)
all_ensg_mart= all_ensg_mart[!is.na(all_ensg_mart\$entrezgene),]
all_ensg_mart = all_ensg_mart[!duplicated(all_ensg_mart\$entrezgene), ]

idx_${i}=match(all_ensg_mart\$ensembl_gene_id, sub(\"\\\\..+\",\"\",rownames(${i}\$table)))
exp.fc=${i}\$table\$logFC[idx_${i}]
names(exp.fc)=all_ensg_mart\$entrezgene[idx_${i}]
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
write.table(rbind(fc.kegg.p\$greater, fc.kegg.p\$less), file = \"${i}_fc.kegg.p.tsv\", sep = \"\\t\")
fc.kegg.sig=sigGeneSet(fc.kegg.p)
write.table(rbind(fc.kegg.sig$greater, fc.kegg.sig$less), file = \"${i}.kegg.sig.tsv\", sep =\"\\t\")

fc.kegg.p.2p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL,same.dir = F)
write.table(rbind(fc.kegg.p.2p\$greater, fc.kegg.p.2p\$less), file = \"${i}_fc.kegg.p.2p.tsv\", sep = \"\\t\")
fc.kegg.2p.sig=sigGeneSet(fc.kegg.p.2p)
write.table(rbind(fc.kegg.2p.sig$greater, fc.kegg.2p.sig$less), file = \"${i}.kegg.2d.sig.tsv\", sep =\"\\t\")


sel <- fc.kegg.p\$greater[, \"q.val\"] < 0.1 & "'!'"is.na(fc.kegg.p\$greater[, \"q.val\"])
path.ids <- rownames(fc.kegg.p\$greater)[sel]
sel.l <- fc.kegg.p\$less[, \"q.val\"] < 0.1 & "'!'"is.na(fc.kegg.p\$less[,\"q.val\"])
path.ids.l <- rownames(fc.kegg.p\$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2[1:20], function(pid) pathview( gene.data = exp.fc, pathway.id = pid, species = \"hsa\", out.suffix=\"${i}\"))


sel <- fc.kegg.p.2p\$greater[, \"q.val\"] < 0.1 & "'!'"is.na(fc.kegg.p.2p\$greater[, \"q.val\"])
path.ids <- rownames(fc.kegg.p.2p\$greater)[sel]
sel.l <- fc.kegg.p.2p\$less[, \"q.val\"] < 0.1 & "'!'"is.na(fc.kegg.p.2p\$less[,\"q.val\"])
path.ids.l <- rownames(fc.kegg.p.2p\$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2[1:20], function(pid) pathview( gene.data = exp.fc, pathway.id = pid, species = \"hsa\", out.suffix=\"2p.${i}\"))

${i}.kegg.esg.up <- esset.grp(${i}.kegg.p$greater, ${i}, gsets = kegg.gs, ref = hn, samp = dcis, test4up = T, output = T, outname = \"${i}.kegg.up\", make.plot = F)
${i}.kegg.esg.dn <- esset.grp(${i}.kegg.p$less, ${i}, gsets = kegg.gs, ref = hn, samp = dcis,test4up = F, output = T, outname = \"${i}.kegg.dn\", make.plot = F)

"; done




###trans
for i in `cat test` ; do echo "

summary(dt_${i} <- decideTestsDGE(${i}, p=0.01, adjust=\"BH\"))
isDE_${i} <- as.logical(dt_${i})
DE_${i} <- rownames(y)[isDE_${i}]
${i}\$table\$FDR <- p.adjust(${i}\$table\$PValue, method=\"BH\")
DEG_${i} <- ${i}\$table[${i}\$table\$FDR < 0.01,]
DEG_${i}=DEG_${i}[DEG_${i}\$FDR > 0,]
DEG_${i}=DEG_${i}[DEG_${i}\$PValue < 0.01,]
#DEG_${i}=DEG_${i}[(DEG_${i}\$logFC<(-2))|(DEG_${i}\$logFC>2),]

results_DEG_${i} <- getBM(attributes = c(\"ensembl_gene_id\", \"ensembl_transcript_id\", \"hgnc_symbol\", \"entrezgene\", \"gene_biotype\", \"description\", \"chromosome_name\", \"start_position\", \"end_position\", \"strand\"), filters = \"ensembl_transcript_id\", values = sub(\"\\\\..+\",\"\",rownames(DEG_${i})),mart = mart)


idx_DEG_${i}=match(sub(\"\\\\..+\",\"\",rownames(DEG_${i})),results_DEG_${i}\$ensembl_transcript_id)

DEG_${i}\$ensembl_gene_id=results_DEG_${i}\$ensembl_gene_id[idx_DEG_${i}]
DEG_${i}\$ensembl_transcript_id=results_DEG_${i}\$ensembl_transcript_id[idx_DEG_${i}]
DEG_${i}\$hgnc_symbol=results_DEG_${i}\$hgnc_symbol[idx_DEG_${i}]
DEG_${i}\$entrezgene=results_DEG_${i}\$entrezgene[idx_DEG_${i}]
DEG_${i}\$gene_biotype=results_DEG_${i}\$gene_biotype[idx_DEG_${i}]
DEG_${i}\$description=results_DEG_${i}\$description[idx_DEG_${i}]
DEG_${i}\$chromosome_name=results_DEG_${i}\$chromosome_name[idx_DEG_${i}]
DEG_${i}\$start_position=results_DEG_${i}\$start_position[idx_DEG_${i}]
DEG_${i}\$end_position=results_DEG_${i}\$end_position[idx_DEG_${i}]
DEG_${i}\$strand=results_DEG_${i}\$strand[idx_DEG_${i}]
####DEG_${i}\$go_id=results_DEG_${i}\$go_id[idx_DEG_${i}]


write.table(DEG_${i}, file=\"${i}_DEG.tsv\", sep=\"\\t\")
top_${i}=rownames(DEG_${i})
${i}_cpm <- cpm(y)[top_${i}, ]
z_score_${i} <- ((${i}_cpm - rowMeans(${i}_cpm))/apply(${i}_cpm,1,sd))

write.table(${i}_cpm, file=\"${i}_DEG_cpm.tsv\", sep=\"\\t\", row.names =T)

unique_DEG_${i}=unique(DEG_${i}\$entrezgene)


pdf(\"${i}_heatmap.pdf\")
pheatmap(z_score_${i}, scale=\"row\", annotation = annotation, annotation_colors = ann_col, cluster_cols=T, cluster_rows=T, show_rownames=F,color = colorRampPalette(rev(brewer.pal(n = 11, name = \"RdYlBu\")))(100))
dev.off()

pdf(\"${i}_MA_plot.pdf\")
plotSmear(${i}, de.tags=DE_${i})
abline(h=c(-1,1), col=\"blue\")
dev.off()

pdf(\"${i}_pvalues.pdf\")
hist(${i}\$table\$PValue, breaks=seq(0, 1, 0.05))
dev.off()

pdf(\"${i}_lrt_RLE.pdf\")
plotRLE(${i}\$fitted.values, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

anno_bp <- annFUN.org(\"BP\", mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
allGenes_bp <- unique(unlist(anno_bp))
uniqueDEG_${i} <- unique(results_DEG_${i}\$ensembl_gene_id)
geneList_${i} <- factor(as.integer(allGenes_bp %in% uniqueDEG_${i}))
names(geneList_${i}) <- allGenes_bp

GOdata_bp_${i} <- new(\"topGOdata\", ontology = \"BP\", allGenes = geneList_${i}, nodeSize = 10, annot = annFUN.org, mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
resultFisher_bp_${i} <- runTest(GOdata_bp_${i}, algorithm = \"classic\", statistic = \"fisher\")
resultKS_bp_${i} <- runTest(GOdata_bp_${i}, algorithm = \"classic\", statistic = \"ks\")
resultKS.elim_bp_${i} <- runTest(GOdata_bp_${i}, algorithm = \"elim\", statistic = \"ks\")

pdf(\"${i}_gotree_BP.pdf\")
showSigOfNodes(GOdata_bp_${i}, score(resultFisher_bp_${i}), firstSigNodes = 10, useInfo = \"all\")
dev.off()

allRes_BP_${i} <- GenTable(GOdata_bp_${i}, classicFisher = resultFisher_bp_${i}, classicKS = resultKS_bp_${i}, elimKS = resultKS.elim_bp_${i}, orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 30)
for (i in 1:nrow(allRes_BP_${i})){ allRes_BP_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(allRes_BP_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(allRes_BP_${i})){ allRes_BP_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(allRes_BP_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(allRes_BP_${i}, file=\"${i}_top_BP.tsv\", sep=\"\\t\", row.names =F)

anno_cc <- annFUN.org(\"CC\", mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
allGenes_cc <- unique(unlist(anno_cc))
uniqueDEG_${i} <- unique(results_DEG_${i}\$ensembl_gene_id)
geneList_${i} <- factor(as.integer(allGenes_cc %in% uniqueDEG_${i}))
names(geneList_${i}) <- allGenes_cc

GOdata_cc_${i} <- new(\"topGOdata\", ontology = \"CC\", allGenes = geneList_${i}, nodeSize = 10, annot = annFUN.org, mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
resultFisher_cc_${i} <- runTest(GOdata_cc_${i}, algorithm = \"classic\", statistic = \"fisher\")
resultKS_cc_${i} <- runTest(GOdata_cc_${i}, algorithm = \"classic\", statistic = \"ks\")
resultKS.elim_cc_${i} <- runTest(GOdata_cc_${i}, algorithm = \"elim\", statistic = \"ks\")

pdf(\"${i}_gotree_CC.pdf\")
showSigOfNodes(GOdata_cc_${i}, score(resultFisher_cc_${i}), firstSigNodes = 10, useInfo = \"all\")
dev.off()

allRes_cc_${i} <- GenTable(GOdata_cc_${i}, classicFisher = resultFisher_cc_${i}, classicKS = resultKS_cc_${i}, elimKS = resultKS.elim_cc_${i}, orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 30)
for (i in 1:nrow(allRes_cc_${i})){ allRes_cc_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(allRes_cc_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(allRes_cc_${i})){ allRes_cc_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(allRes_cc_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(allRes_cc_${i}, file=\"${i}_top_CC.tsv\", sep=\"\\t\", row.names =F)

anno_mf <- annFUN.org(\"MF\", mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
allGenes_mf <- unique(unlist(anno_mf))
uniqueDEG_${i} <- unique(results_DEG_${i}\$ensembl_gene_id)
geneList_${i} <- factor(as.integer(allGenes_mf %in% uniqueDEG_${i}))
names(geneList_${i}) <- allGenes_mf

GOdata_mf_${i} <- new(\"topGOdata\", ontology = \"MF\", allGenes = geneList_${i}, nodeSize = 10, annot = annFUN.org, mapping = \"org.Hs.eg.db\", ID = \"ensembl\")
resultFisher_mf_${i} <- runTest(GOdata_mf_${i}, algorithm = \"classic\", statistic = \"fisher\")
resultKS_mf_${i} <- runTest(GOdata_mf_${i}, algorithm = \"classic\", statistic = \"ks\")
resultKS.elim_mf_${i} <- runTest(GOdata_mf_${i}, algorithm = \"elim\", statistic = \"ks\")

pdf(\"${i}_gotree_MF.pdf\")
showSigOfNodes(GOdata_mf_${i}, score(resultFisher_mf_${i}), firstSigNodes = 10, useInfo = \"all\")
dev.off()

allRes_mf_${i} <- GenTable(GOdata_mf_${i}, classicFisher = resultFisher_mf_${i}, classicKS = resultKS_mf_${i}, elimKS = resultKS.elim_mf_${i}, orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 30)
for (i in 1:nrow(allRes_mf_${i})){ allRes_mf_${i}\$hgnc_symbol[i]=paste(as.vector(results_DEG_${i}\$hgnc_symbol)[which(grepl(allRes_mf_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
for (i in 1:nrow(allRes_mf_${i})){ allRes_mf_${i}\$ensg[i]=paste(as.vector(results_DEG_${i}\$ensembl_gene_id)[which(grepl(allRes_mf_${i}\$GO.ID[i],results_DEG_${i}\$go_id))],collapse=\",\")}
write.table(allRes_mf_${i}, file=\"${i}_top_MF.tsv\", sep=\"\\t\", row.names =F)

unique_DEG_${i}=unique(DEG_${i}\$entrezgene)
go_DEG_${i}=goana(unique_DEG_${i})
top_BP_DEG_${i}=topGO(go_DEG_${i},ont=\"BP\",n=30)
top_MF_DEG_${i}=topGO(go_DEG_${i},ont=\"MF\",n=30)
top_CC_DEG_${i}=topGO(go_DEG_${i},ont=\"CC\",n=30)

write.table(top_BP_DEG_${i}, file=\"${i}_top_BP.tsv\", sep=\"\\t\")
write.table(top_MF_DEG_${i}, file=\"${i}_top_MF.tsv\", sep=\"\\t\")
write.table(top_CC_DEG_${i}, file=\"${i}_top_CC.tsv\", sep=\"\\t\")

"; done

