#Comp Project

library(tximport)
library(readr)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

#load in outputs
samples_df<- data.frame(
  sample_name = c("GRN 24","GRN 66","GRN 31","GRN 2","GRN 28","GRN 41", "Control 30", "Control 12", "Control 64", "Control 1","C9 39","C9 40","C9 33","C9 11","C9 27","C9 42", "MAPT 43", "MAPT 58", "MAPT 59", "MAPT 62", "MAPT 34","MAPT 18"),
  patient_ID = c("Rimod_24","Rimod_66","Rimod_31","Rimod_2","Rimod_28","Rimod_41", "Rimod_30", "Rimod_12", "Rimod_64", "Rimod_1","Rimod_39","Rimod_40","Rimod_33","Rimod_11","Rimod_27","Rimod_42", "Rimod_43", "Rimod_58", "Rimod_59", "Rimod_62", "Rimod_34","Rimod_18"),
  condition = c(rep("GRN Mutation",6), rep("Control",4), rep("C9orf72",6),rep("MAPT Mutation",6)), 
  group_folder = c(rep("GRN",6),rep("Controls",4), rep("C9 mutation", 6),rep("MAPT",  6))
)


#IMPORTANT: Update filepath prior to running
samples_df$filepath<-file.path(
  "C:/Users/Lauren/Documents/Comp_Project",
  samples_df$group_folder,
  samples_df$patient_ID,
  "output_folder",
  "abundance.tsv"
)

files<- samples_df$filepath
names(files) <- samples_df$sample_name

#gene annotation
#IMPORTANT: Update filepath before running
library(GenomicFeatures)
gtf_path <- "C:/Users/Lauren/Documents/Comp_Project/gencode.v47.annotation.gtf"

txdb <- makeTxDbFromGFF(gtf_path)

tx2gene <- select(txdb, keys = keys(txdb, "TXNAME"), columns = c("TXNAME", "GENEID"), keytype = "TXNAME")
head(tx2gene)


txi<- tximport(files = files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples_df, design = ~condition)

#map ensembl ids
rownames(dds)<-sub("\\..*","",rownames(dds))

gene_symbols<-mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL",keytype = "ENSEMBL", multiVals = "first")

rowData(dds)$symbol<-gene_symbols

#sort mitochondria genes
mito_genes<- grep("^MT-", gene_symbols,value=TRUE)

dds_mito<-dds[which(rowData(dds)$symbol %in% mito_genes),]

#set control as reference
dds_mito$condition<- relevel(dds_mito$condition, ref = "Control")
dds_mito<-DESeq(dds_mito)

#results by mutation
res_C9_mito<-results(dds_mito, contrast =c("condition", "C9orf72","Control"))
head(res_C9_mito)
summary(res_C9_mito)


res_GRN_mito<-results(dds_mito, contrast =c("condition", "GRN Mutation","Control"))
head(res_GRN_mito)
summary(res_GRN_mito)

res_MAPT_mito<-results(dds_mito, contrast =c("condition", "MAPT Mutation","Control"))
head(res_MAPT_mito)
summary(res_MAPT_mito)

#filter for significant genes
sig_mito_C9<-res_C9_mito[which(res_C9_mito$padj<0.1),]
sig_mito_C9

sig_mito_GRN<-res_GRN_mito[which(res_GRN_mito$padj<0.1),]
sig_mito_GRN

sig_mito_MAPT<-res_MAPT_mito[which(res_MAPT_mito$padj<0.1),]
sig_mito_MAPT


#plots:
library(ggplot2)


#volcano plots
res_C9_mito$significant <- ifelse(res_C9_mito$padj<0.1,"Significant","Not Significant")
res_GRN_mito$significant <- ifelse(res_GRN_mito$padj<0.1,"Significant","Not Significant")
res_MAPT_mito$significant <- ifelse(res_MAPT_mito$padj<0.1,"Significant","Not Significant")

ggplot(as.data.frame(res_C9_mito), aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color = significant))+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+
  labs(
    title = "C9orf72 vs Control (Mitochondrial Genes)",
    x="Log2 Fold Change",
    y= "-log10(p-value)"
  )
ggplot(as.data.frame(res_GRN_mito), aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color = significant))+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+
  labs(
    title = "GRN vs Control (Mitochondrial Genes)",
    x="Log2 Fold Change",
    y= "-log10(p-value)"
  )

ggplot(as.data.frame(res_MAPT_mito), aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color = significant))+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+
  labs(
    title = "MAPT vs Control (Mitochondrial Genes)",
    x="Log2 Fold Change",
    y= "-log10(p-value)"
  )


#boxplots:
padj.cutoff<-0.1
p.cutoff<- 0.1

C9up<-sum(res_C9_mito$pvalue < p.cutoff & res_C9_mito$log2FoldChange> 0, na.rm=TRUE)
C9down<-sum(res_C9_mito$pvalue < p.cutoff & res_C9_mito$log2FoldChange< 0, na.rm=TRUE)

GRNup<-sum(res_GRN_mito$pvalue < p.cutoff & res_GRN_mito$log2FoldChange> 0, na.rm=TRUE)
GRNdown<-sum(res_GRN_mito$pvalue < p.cutoff & res_GRN_mito$log2FoldChange< 0, na.rm=TRUE)

MAPTup<-sum(res_MAPT_mito$padj < p.cutoff & res_MAPT_mito$log2FoldChange> 0, na.rm=TRUE)
MAPTdown<-sum(res_MAPT_mito$padj < p.cutoff & res_MAPT_mito$log2FoldChange< 0, na.rm=TRUE)

df_counts <- data.frame(
  group = rep(c("C9orf72", "GRN Mutation", "MAPT Mutation"), each = 2),
  direction = rep(c("Up", "Down"), times=3),
  count = c(C9up, C9down, GRNup, GRNdown, MAPTup,MAPTdown)
)

ggplot(df_counts, aes(x=group, y = count, fill=direction))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(
    name="Direction",
    values = c("Up"="navy",
               "Down"="darkred")
  )+
  labs(
    title = "Significant Mitochondrial Gene Dysregulation",
    x="",
    y= "Number of genes (padj<0.1)"
  )+
  theme_minimal(base_size = 14)+
  theme(
    axis.text.x = element_text(face="bold"),
    legend.position = "top"
  )


#fold change dist
ggplot(as.data.frame(res_C9_mito),aes(x=log2FoldChange))+
  geom_density()+
  labs(title="C9orf72 Mitochondrial FC Distribution")

ggplot(as.data.frame(res_GRN_mito),aes(x=log2FoldChange))+
  geom_density()+
  labs(title="GRN Mitochondrial FC Distribution")

ggplot(as.data.frame(res_MAPT_mito),aes(x=log2FoldChange))+
  geom_density()+
  labs(title="MAPT Mitochondrial FC Distribution")


#heatmaps
library(pheatmap)

vsd_mito <- varianceStabilizingTransformation(dds_mito, blind = TRUE)

annotation <- data.frame(condition = samples_df$condition)
rownames(annotation) <- samples_df$sample_name

all_C9_ids  <- rownames(res_C9_mito)
all_GRN_ids <- rownames(res_GRN_mito)
all_MAPT_ids <- rownames(res_MAPT_mito)

ranked_C9  <- all_C9_ids[order(-abs(res_C9_mito$log2FoldChange))]
ranked_GRN <- all_GRN_ids[order(-abs(res_GRN_mito$log2FoldChange))]
ranked_MAPT <-all_MAPT_ids[order(-abs(res_MAPT_mito$log2FoldChange))]

plot_C9_ids  <- intersect(head(ranked_C9, 20),  rownames(vsd_mito))
plot_GRN_ids <- intersect(head(ranked_GRN,20),   rownames(vsd_mito))
plot_MAPT_ids <- intersect(head(ranked_MAPT,20),   rownames(vsd_mito))

# stop if broken
if(length(plot_C9_ids)==0) stop("No C9 IDs found in VST object")
if(length(plot_GRN_ids)==0) stop("No GRN IDs found in VST object")
if(length(plot_MAPT_ids)==0) stop("No MAPT IDs found in VST object")

mat_C9  <- assay(vsd_mito)[plot_C9_ids,  colnames(vsd_mito)]
mat_GRN <- assay(vsd_mito)[plot_GRN_ids, colnames(vsd_mito)]
mat_MAPT <- assay(vsd_mito)[plot_MAPT_ids, colnames(vsd_mito)]

pheatmap(mat_C9,
         annotation_col = annotation,
         scale          = "row",
         main           = "Top 20 Mito Genes: C9orf72 vs Control")

pheatmap(mat_GRN,
         annotation_col = annotation,
         scale          = "row",
         main           = "Top 20 Mito Genes: GRN Mutation vs Control")

pheatmap(mat_MAPT,
         annotation_col = annotation,
         scale          = "row",
         main           = "Top 20 Mito Genes: MAPT Mutation vs Control")


#PCA Plot
plotPCA(vsd_mito, intgroup="condition") + ggtitle("PCA: Mito Genes")



