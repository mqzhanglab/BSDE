# Introduction

This markdown is for individual case-control study based on Barycenter and Wasserstein distance estimation on the gene differential expression (DE) analysis of certain cell type based on the idiopathic pulmonary fibrosis.

# Initial Setup

R packages required in this pipeline:
```{r load_libraries,eval=runcode}
library("MASS")
library("emdbook")
library("moments")
library("MAST")
library("lme4")
library("DESeq2")
library("abind")
library("doRNG")
library("reticulate")
library("BSDE")


```

# Computation

Here we present a demo version of our real data analysis,We take the basic parameters, i.e., the mean, dispersion and dropout parameter from the distribution of the real reference scRNAseq database from paper [Single-cell RNA sequencing identifies diverse roles of epithelial cells in idiopathic pulmonary fibrosis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5135277/).

We focus on the cell type SCGB3A2+ and limit our analysis in the first 1000 genes.The trimmed file are at the data/IPF folder.

## Basic Settings

Senario settings in our analysis are as follows: 


```{r basic_settings,eval=runcode}
#senarios
covariate_flag=NA #c(NA, "quantile99")
perm_label=0 #perm_label =0 means calculate the observed data other wise, permutated data
merge_method="python" # c("R","python","python_exp")
lambda_fold=2


dir = "~/Dropbox/Barycenter/BSDEpackage/BSDE/data/"
setwd(dir)
```

## Data Preparation
```{R data_load,eval=runcode}
###########input###############
meta=readRDS(paste0(dir,"IPF/meta.rds"))

#input phenotype
cur_individual=unique(meta$individual)

#input bulk
rawM_bulk=readRDS(paste0(dir,"IPF/matrix_bulk.rds"))
rawM_cell=readRDS(paste0(dir,"IPF/ind_cellnum.rds"))
rawM_bulk=t(apply(rawM_bulk,1,function(x) x*100/rawM_cell)) 
colnames(rawM_bulk)=names(rawM_cell)
logcount_data_bulk=log2(rawM_bulk+1)

#input counts
rawM=readRDS(paste0(dir,"IPF/rawM.rds"))
logcount_data=log2(rawM+1)
dim(logcount_data)
dim(meta)
meta_ind=meta[match(unique(meta$individual),meta$individual),c("individual","diagnosis")]

rdp=readRDS(paste0(dir,"IPF/ind_readdepth.rds"))

meta_ind$rdp=scale(rdp)
meta_ind$individual=scale(as.numeric(meta_ind$individual))
meta_ind$diagnosis=meta_ind$diagnosis


```

After read in the data, we apply inverse probability weight to remove the covariance effects. 
Here we use the factor rdp as the ipw denominator example, in realworld, the denominator should be individual related factors such as age and sex, rather than experiment related ones.

```{r ipw,eval=runcode,message = FALSE,warning = FALSE,include = FALSE,results = FALSE}
library("ipw")

temp=ipwpoint(exposure=diagnosis,family="binomial", link="logit",numerator=~1,denominator = ~ rdp,data=meta_ind)
plot(temp$ipw.weights, rdp)
ipww=temp$ipw.weights
ipww

```


Then we run our method for that.
```{r run_BSDE,eval=FALSE}

dim(logcount_data)
logcount_data[1:10,1:10]
cell_id=meta$cell
gene_id=dimnames(logcount_data)[[1]]
rownames(logcount_data)=gene_id
colnames(logcount_data)=cell_id

op_pvalw=matrix(ncol=1,nrow=nrow(logcount_data))

bc_case_ob=list()
bc_ctrl_ob=list()
bcw_case_ob=list()
bcw_ctrl_ob=list()
#for(ig in 1:nrow(logcount_data)){
for(ig in 1:10){
  
  bc_case_ob[[ig]]=list()
  bc_ctrl_ob[[ig]]=list()
  bcw_case_ob[[ig]]=list()
  bcw_ctrl_ob[[ig]]=list()
  
  cur_sim=logcount_data[ig,]
  cur_ind=meta$individual
  cur_pheno=meta$diagnosis
  
  tempa=tryCatch({cal_w2_pval(count_per_gene=cur_sim,meta_individual=cur_ind,meta_phenotype=cur_pheno,perm_num=200,unif_round_unit=0.5,l_fold=lambda_fold,merge_method=merge_method,weight=ipww,shrink=FALSE)}, error = function(e) {NA} )
  op_pvalw[ig]=tryCatch({tempa[[1]]}, error = function(e) {NA} )
  bcw_case_ob[[ig]]=tryCatch({tempa[[2]]}, error = function(e) {NA} )
  bcw_ctrl_ob[[ig]]=tryCatch({tempa[[3]]}, error = function(e) {NA} )
  if(!(is.na(op_pvalw[ig])) && op_pvalw[ig]<2/200){
    op_pvalw[ig]=tryCatch({cal_w2_pval(count_per_gene=cur_sim,meta_individual=cur_ind,meta_phenotype=cur_pheno,perm_num=1000,unif_round_unit=0.5,l_fold=lambda_fold,merge_method=merge_method,weight=ipww,shrink=FALSE)[[1]]}, error = function(e) {NA} )}
  if(!(is.na(op_pvalw[ig])) && op_pvalw[ig]<2/1000){
    op_pvalw[ig]=tryCatch({cal_w2_pval(count_per_gene=cur_sim,meta_individual=cur_ind,meta_phenotype=cur_pheno,perm_num=5000,unif_round_unit=0.5,l_fold=lambda_fold,merge_method=merge_method,weight=ipww,shrink=FALSE)[[1]]}, error = function(e) {NA} )}
}

names(op_pvalw)=gene_id


```

# Result Visualization

Once got the p-values, we can present our results as follows.

## Data Preparation

We've done our the analysis based on the top 3000 expressed genes among all clusters, and use FDR method multi-effect correction. Then we will do some data analysis and visualization.


```{r plot_functions,eval=runcode}
n_cluster=31
case_col="#00539B"
ctrl_col="#a3d6f0" #"#4B9CD3"

set.seed(1)
library(ggplot2)
library(stringr)
library(scales)

#geom_split_violin function
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

t_col = function(color, percent = 60, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
```

```{r result_load,eval=runcode}
#read_in
opw_padj=readRDS(paste0(dir,"IPF/opw_padj.rds"))
gene_id=names(opw_padj[[1]])
tmeta=readRDS(paste0(dir,"IPF/tmeta.rds"))
cell_num_table=readRDS(paste0(dir,"IPF/cell_num_table.rds"))
cluster=unique(tmeta$cluster)
total_individual=unique(tmeta$individual)

cur_cluster="Proliferating Epithelial Cells" #Proliferating Epithelial Cells  #Plasma Cells #Endothelial Cells
cluster=unique(tmeta$celltype)
cluster_plot=gsub("Proliferating","P.",cluster)
cluster_plot=gsub("Differentiating","Diff.",cluster_plot)
cluster_plot=gsub("Transitional","Trans.",cluster_plot)
cluster_plot=gsub("+ SCGB1A","/",cluster_plot)

```


## Boxplot of the number of cells per subjects.

```{r cell_per_subject,eval=runcode}
pval_num=sapply(opw_padj,function(x){sum(x<=0.1,na.rm=TRUE)}) #adjusted FDR <0.1

population=tmeta$population[match(cluster,tmeta$cluster)]
sub_cluster_index=which(population=="Epithelial")

diagnose=tmeta$Status[match(total_individual,tmeta$individual)]
cell_num_table=cell_num_table[sub_cluster_index,]
pval_num=pval_num[sub_cluster_index]
porder=order(pval_num,decreasing = TRUE)
this_table=data.frame(cell_num=as.numeric(cell_num_table),
                      cluster=factor(rep(cluster_plot[sub_cluster_index],times=ncol(cell_num_table))
                      ,levels=cluster_plot[sub_cluster_index][porder]),
                      Diagnose=factor(rep(diagnose,each=nrow(cell_num_table)),levels=c("ILD","Control")))
this_table$cur_col=rep(case_col,nrow(this_table))
this_table$cur_col[this_table$diagnose=="Control"]=ctrl_col

this_table=this_table[this_table$cell_num>0,]
#violin plot cell number
cols = hue_pal()(length(unique(this_table$cluster)))
cols2=sapply(cols,t_col)

p1=ggplot(this_table, aes(x=cluster, y=cell_num, fill=Diagnose)) +
  #scale_fill_manual(values=hue_pal()(2)) +
  scale_fill_manual(values=c(case_col,ctrl_col)) +
  labs(title="Number of Epithelial Cells per Subject",x=" ", y = "Number of Cells", tag = " ")+
  #geom_split_violin(trim=FALSE)+#geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot()+ scale_y_continuous(trans='log10') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p1)

```

## bar plot of number of differentially expressed genes per subjects.
```{r pval_enriched,eval=runcode}
this_table=data.frame(pval_num=pval_num,
                      cluster=factor(cluster_plot[sub_cluster_index],
                      levels=cluster_plot[sub_cluster_index][porder]))

p2=ggplot(this_table, aes(x=cluster, y=pval_num,fill=1)) +
  labs(title="Number of DEGs in Epithelial Cells",x=" ", y = "Number of DEGs", tag = " ")+
  geom_bar(stat="identity")+ #geom_jitter(shape=16, position=position_jitter(0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p2)
```

## Gene set enrichment analysis

We also perform gene set enrichment analysis using the [GOSeq package](https://bioconductor.org/packages/release/bioc/html/goseq.html).For manually annotation, we download gencode.v37.annotation.gtf from [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz).


```{r goseq,eval=runcode,message = FALSE,warning = FALSE,include = FALSE,results = FALSE}
#go analysis
library(goseq)
require("biomaRt")
library("biomaRt")
library("org.Hs.eg.db") # remember to install it if you don't have it already
library("goseq")
library(GO.db)
# We also need to extract other information from gtf files since they are not included in txdb:
gencode_gtf = rtracklayer::import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz")
gencode_gtf = gencode_gtf[!duplicated(gencode_gtf$gene_id), ]

go_table=matrix(nrow=0,ncol=3)
colnames(go_table)=c("cluster","term","nlog10pval")

for(i_cluster in sub_cluster_index[porder]){
  
  #KEGG analysis
  cur_pval=opw_padj[[i_cluster]]
  genes=as.integer(cur_pval<=0.05)
  names(genes)=names(cur_pval)
  
  #match
  gencode_index=match(names(cur_pval), gencode_gtf$gene_name)
  #gencode_gtf = gencode_gtf[gencode_index[!is.na(gencode_index)],
  #                          c("gene_id", "source", "gene_type", "gene_name")]
  
  genes=genes[!is.na(gencode_index)]
  genes=genes[!is.na(genes)]
  #unique the names
  gene1=unique(names(genes[genes==1]))
  genesa=unique(names(genes))
  genes=rep(0,length(genesa))
  names(genes)=genesa
  genes[match(gene1,genesa)]=1
  # gene_name is gene common name, for example, BRCA
  # gene_id is gene ensembl id, ENSG0000000001111
  
  #GO
  pwf=nullp(genes,'hg19','geneSymbol')
  GO.wall=goseq(pwf,"hg19","geneSymbol")
  #head(GO.wall)
  sort(GO.wall$over_represented_pvalue)[1:10]
  
  #multiple correction
  GO.wall$qvalue=p.adjust(GO.wall$over_represented_pvalue,method="fdr")
  enriched.GO = GO.wall$category[GO.wall$qvalue<.2]
  
  #print out the significant GOs
  n_sig=length(enriched.GO)
  if(n_sig>0){
    #print(GO.wall[1:min(5,n_sig),c(1,2,6,7)])
    this_table=data.frame(cluster=rep(cluster[i_cluster],n_sig),
                          term=GO.wall$term[1:n_sig],
                          nlog10pval=-log10(GO.wall$qvalue[1:n_sig]))
    this_table=this_table[1:min(n_sig,5),]
    go_table=rbind(go_table,this_table)
  }
}
go_table$label=apply(go_table[,c("term","cluster")],1,
                     function(x){paste0(x[1],"    ", str_pad(x[2],width=11,side="left"))})
go_table$label=factor(go_table$label,levels=rev(go_table$label))
```

```{r goseq_plot,eval=runcode}
p3=ggplot(go_table, aes(x=label, y=nlog10pval,fill=cluster)) +
  labs(title="GO: DEGs Enriched Pathways",x=" ", y = "-log10 Enrichment q-values", tag = " ")+
  geom_bar(stat="identity")+ #geom_jitter(shape=16, position=position_jitter(0.2))+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p3)

```


## Illustration of the differentially expressed genes

```{r sig_genes,eval=runcode}

op=par(mfrow=c(5,4),oma=c(0,0,0,0),mar=c(3,3,3,0),cex.lab=2,cex.axis=2,cex.main=2)
#focus on "SCGB3A2+","AT2","Basal","MUC5B+","AT1"
for(i_cluster in c(7,11,2,8,12)){ 
  count=0 #initiate
  sig_names=matrix(ncol=1,nrow=0)
  
  cur_cluster=as.character(unique(tmeta$cluster)[i_cluster])
  meta=tmeta[tmeta$cluster==cur_cluster,]
  meta$diagnosis=as.numeric(meta$diagnosis=="Control")
  cur_individual=unique(meta$individual)
  
  cur_cluster=cluster[i_cluster]
  cur_pval=NA
  cur_padj=NA
  
  cur_pval=readRDS(paste0(dir,"IPF/op_pval/op_pvalw.rds"))[[i_cluster]]
  names(cur_pval)=gene_id
  names(cur_pval)=gsub("ABR_ENSG00000159842","ABR",names(cur_pval))
  cur_padj=p.adjust(cur_pval,method="fdr")
  
  bcw_case=readRDS(paste0(dir,"IPF/op_pval/bcw_case.rds"))[[i_cluster]]
  bcw_ctrl=readRDS(paste0(dir,"IPF/op_pval/bcw_ctrl.rds"))[[i_cluster]]
  
  cur_n=100
  
  sig_index_opw=c(order(cur_pval,decreasing = FALSE))[1:sum(cur_padj<=0.1,na.rm = TRUE)]
  
  mt_index=grep("^MT",names(cur_pval)) #remove MT genes
  sig_index_opw=setdiff(sig_index_opw,mt_index)
  
  ###plot 
  for(i_sig in c(sig_index_opw)){
    
    cur_case=bcw_case[[i_sig]]
    cur_ctrl=bcw_ctrl[[i_sig]]
    
    case_dx=as.numeric(names(cur_case))
    case_dy=as.numeric(cur_case)
    ctrl_dx=as.numeric(names(cur_ctrl))
    ctrl_dy=as.numeric(cur_ctrl)
    
    if(is.na(match(names(cur_pval)[i_sig],sig_names)) &&
       (sum(case_dy[case_dx>=1]) + sum(ctrl_dy[ctrl_dx>=1])) >=0.5){ #total proportion >=1 more than 40%
      
      plot(case_dx, case_dy, type = "l",xlab="",ylab="",
           xlim=c(0,max(case_dx,ctrl_dx)),ylim=c(0,max(case_dy,ctrl_dy)),
           main=gene_id[i_sig],frame.plot = FALSE, axes.y=FALSE)#, axes=FALSE)
      lines(ctrl_dx, ctrl_dy, type = "l")
      polygon(c(case_dx,rev(case_dx)), c(case_dy,rep(0,length(case_dy))), col=t_col(case_col,percent = 20), border=NA)
      polygon(c(ctrl_dx,rev(ctrl_dx)), c(ctrl_dy,rep(0,length(ctrl_dy))), col=t_col(ctrl_col,percent = 50), border=NA)
      
      count=count+1
      sig_names=c(sig_names,names(cur_pval)[i_sig])
      #print(c(names(cur_pval)[i_sig],cluster[i_cluster],count))
    }
    if(count==4){break}
  }
}
par(op)

```

