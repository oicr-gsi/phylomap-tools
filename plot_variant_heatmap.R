#
# Plot a phylogenetic tree with associated mutations
#

#
# this method of aligning the tree and tiles is from:
# https://thackl.github.io/ggtree-composite-plots
#



require(optparse)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default="results", 
              help="data directory", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default="tree.nwk", 
              help="newick tree file", metavar="character"),
  make_option(c("-v", "--variants"), type="character", default="alleles.tsv", 
              help="variants file", metavar="character"),
  make_option(c("-c", "--covariates"), type="character", default=NULL, 
              help="covariate file", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="default", 
              help="output prefix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="plots", 
              help="output file name [default= %default]", metavar="character"),
  
  ### this option is very specific to the metadata, can i specify this in another way, ie. give both key and value for subsetting
  make_option(c("-s", "--sort"), type="character", default=NULL, 
              help="sort order of the heatmap", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#if (is.null(opt$file)){
#  print_help(opt_parser)
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}

#### FOR THIS ANALYSIS
useTree=FALSE
#opt$covariates="covariates.txt"
#opt$prefix="TWHO_heatmap"



require(ggtree)
require(ggplot2)
require(patchwork)
require(reshape2)



# fix the alignment of the tree panel
scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}

### theme for the covariate plots
themes<-list()
themes$covariates<-theme(axis.line = element_blank(), 
                       axis.title.y = element_blank(), 
                       #axis.title.x = element_blank(),
                       axis.title.x = element_text(angle = 90,hjust=0.5),
                       axis.ticks.y = element_blank(), 
                       axis.ticks.x = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       # panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       legend.position = "top",
                       plot.margin=unit(c(0,0,0,0),"cm"))


themes$variants<-theme(axis.line = element_blank(), 
                       axis.title.y = element_blank(), 
                       axis.title.x = element_blank(),
                       axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                       axis.text.x = element_text(angle = 90, hjust = 1),
                       legend.position = "top")

plot_tree<-function(t)
{
  g <- ggtree(t) + 
    #geom_tiplab(align=TRUE) + 
    #xlim(0, 0.00050) +
    scale_y_tree() +
    theme_tree2()
 return(g)
}

plot_variants<-function(a)
{
  a$pos = factor(alleles$pos)
  bases = c("A", "C", "G", "T", "N", "X")
  a$alt_allele = factor(a$alt_allele, bases)
  ambiguous_bases = !(a$alt_allele %in% bases)
  a$alt_allele[ambiguous_bases] = "X"

  cols <- c("blue", "red", "green3", "purple3", "lightgrey", "black")
  g <- ggplot(a, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    theme_bw() + 
    themes$variants + 
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  
  return(g) 
}

plot_variants.2<-function(a)
{
  a$pos = factor(alleles$pos)
  bases = c("A", "C", "G", "T", "N")
  a$alt_allele = factor(a$alt_allele, bases)
  ambiguous_bases = !(a$alt_allele %in% bases)
  a$alt_allele[ambiguous_bases] = "X"
  
  cols <- c("blue", "red", "green3", "purple3", "lightgrey")
  g <- ggplot(a, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    theme_bw() + 
    themes$variants + 
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  
  return(g) 
}



plot_covariate<-function(c,i)
{
  g<- ggplot(c,aes(x=1,y=name)) + 
    geom_tile(aes_string(fill=i)) + 
    xlab(i) +
    themes$covariate +
    scale_fill_discrete(na.translate=FALSE)
  return(g)
}

plot_covariate_text<-function(c,i)
{
  g<- ggplot(c,aes(x=1,y=name)) + 
    #geom_tile(aes(fill="white",colour="white")) + 
    geom_tile(colour="white",fill="white") + 
    geom_text(aes_string(label=i),size=2) +
    xlab(i) +
    themes$covariate +
    theme(legend.position = "none") +
    scale_fill_discrete(na.translate=FALSE)

  return(g)
}

#### MAIN





data_dir<-opt$directory
tree_fn<-opt$tree
covariates_fn<-opt$covariates
alleles_fn<-opt$variants
prefix<-opt$prefix
plot_dir<-opt$out


#### get the allele plot, a heatmap
alleles_path = paste(data_dir, alleles_fn, sep="/")
alleles <- read.table(alleles_path, header=T)
names<-unique(alleles$name)

## covariates
covariates_path<-paste(data_dir,covariates_fn,sep="/")
covariates<-read.table(covariates_path, header=T,as.is=T)
covariate.ids<-names(covariates)[-1] 

### only use covariates that are in the allele table. The allele table will only call if passing completeness threeshold
covariates<-covariates[covariates$name %in% names,]

if(! is.null(opt$site)){
  covariates<-covariates[covariates$Site==opt$site,]
  names<-unique(covariates$name)
  alleles<-alleles[alleles$name %in% names,]
}



dforder<-data.frame(
  name=as.vector(covariates$name),
  pid=as.numeric(matrix(unlist(strsplit(as.vector(covariates$name),".",fixed=T)),ncol=2,byrow=TRUE)[,1]),
  sid=matrix(unlist(strsplit(as.vector(covariates$name),".",fixed=T)),ncol=2,byrow=TRUE)[,2],
  CollectionDate=covariates$CollectionDate,
  Site=covariates$Site,
  Type=covariates$Type,
  OutbreakCase=covariates$OutbreakCase,
  stringsAsFactors=FALSE
  
)
### sort is id, then collection data
dforder<-dforder[order(dforder$pid,dforder$CollectionDate),]

## sort by Type, then Site
dforder<-dforder[order(dforder$Type,dforder$Site,dforder$CollectionDate,dforder$pid),]





### modify the factors so the same order as dforder
alleles$name<-factor(alleles$name, levels=dforder$name)
covariates$name = factor(covariates$name, levels=dforder$name)

### a list of plot panels, ordered from left to right
panels<-list()





### this script SKIPS the tree all together

if(useTree==TRUE){
  ## incorporate the tree, currently done in the plot_variant_tree.R script
}else{
  ## do nothing
  ## tip.order is the current order
  tip.order = dforder$name
}


### get the covariate plots, each a heatmap with no margins
if(! is.null(covariates_fn))
{
  #refdf<-data.frame(name="MN908947.3",location=NA,site=NA,batch=NA)
  
  covariates_path<-paste(data_dir,covariates_fn,sep="/")
  covariates<-read.table(covariates_path, header=T,as.is=T)
  covariate.ids<-names(covariates)[-1] 

  ### only use covariates that are in tip.order
  covariates<-covariates[covariates$name %in% tip.order,]
  ### add covariates that are missing from tip.order
  missing<-tip.order[! tip.order %in% covariates$name]
  #missing.df<-covariates[0,]
  missing.df<-data.frame(matrix(NA,ncol=ncol(covariates),nrow=length(missing)))
  colnames(missing.df)<-colnames(covariates)
  missing.df$name<-missing
  covariates<-rbind(covariates,missing.df)
  
  covariates$name = factor(covariates$name, levels=tip.order)
  
  for (id in covariate.ids){
    if(id=="CollectionDate"){
      panels[[id]]<-plot_covariate_text(covariates,id)
    }else{
      panels[[id]]<-plot_covariate(covariates,id)
    }
  }
}else{
  covariate.ids<-NULL
}


#### get the allele plot, a heatmap
alleles_path = paste(data_dir, alleles_fn, sep="/")
alleles <- read.table(alleles_path, header=T)
alleles$name = factor(alleles$name, levels=tip.order)
#panels$variants<-plot_variants(alleles)   
panels$variants<-plot_variants.2(alleles) 

### set the panel width, 1 for the tree, 0.1 for each covariate and 3 for the variants
if(useTree==TRUE){
  panel.widths<-c(1,rep(0.1,length(covariate.ids)),3)
}else{
  panel.widths<-c(rep(0.1,length(covariate.ids)),3)
}

p<-wrap_plots(panels) + plot_layout(widths = panel.widths, guides="collect")

### data for export
make_table<-function(a,c){
  m<-dcast(a,name~pos,value.var="alt_allele")
  m$name<-as.character(m$name)  ### remove factor
  m[is.na(m)]<-""
  
  ref<-unique(a[,c("pos","ref_allele")])
  ref<-ref[order(ref$pos),]
  refrow<-c("MN908947.3",as.character(ref$ref_allele))
  m<-rbind(refrow,m)
  rownames(m)<-m$name
  m<-m[rev(tip.order),]
  
  rownames(c)<-c$name
  c<-c[rev(tip.order),]
  c<-c[,-1]
  c[is.na(c)]<-""
  mm<-cbind(c,m)
  mm
}
export.table<-make_table(alleles,covariates)
fn<-paste(prefix,"_variant_heatmap.tsv",sep="")
write.table(export.table,file=fn,quote=F,sep="\t",row.names=F,col.names=T)


### save RData
fn<-paste(prefix,"_variant_heatmap.Rdata",sep="")
save.image(file=fn)


### some custom work on this

## save the plot    
#plot_path = sprintf("plots/%s_variant_tree.pdf", prefix)
plot_path = paste(plot_dir,"/",prefix,"_variant_heatmap.pdf",sep="")

# count number of samples, for scaling the plot
num_samples = length(tip.order)
pdf_height = 0.125 * num_samples
pdf_height = max(8, pdf_height)
ggsave(plot_path, p, height=pdf_height, width=20)
