#
# Plot a phylogenetic tree with associated mutations
#

#
# this method of aligning the tree and tiles is from:
# https://thackl.github.io/ggtree-composite-plots
#



require(optparse)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default="results",    help="data directory",          metavar="character"),
  make_option(c("-t", "--tree"),      type="character", default=NULL,         help="newick tree file",        metavar="character"),
  make_option(c("-s", "--sort"),      type="character", default=NULL,         help="sample sort order, by covariate", metavar="character"),
  make_option(c("-v", "--variants"),  type="character", default="alleles.tsv",help="variants file",           metavar="character"),
  make_option(c("-m", "--metadata"),  type="character", default=NULL,         help="covariate metadata file", metavar="character"),
  make_option(c("-c", "--covariates"),type="character", default=NULL,         help="comma separated list of covariates and order to show",metavar="character"),
  make_option("--textcovariates",type="character", default=NULL,       help="comma separated list of covariates that should be displayed as text",metavar="character"),
  make_option(c("-p", "--prefix"),    type="character", default="default",    help="output prefix",           metavar="character"),
  make_option(c("-o", "--out"),       type="character", default="plots",      help="output file name [default= %default]", metavar="character"),
  make_option("--reference_name", type="character", default="MN908947.3", help="the name of the reference", metavar="character"),
  make_option("--title",type="character",default="Genomic Variants",help="a title for the plot",metavar="character"),
  make_option("--height",type="character",default=NULL,help="height of the plot",metavar="character"),
  make_option("--width",type="character",default=NULL,help="width of the plot",metavar="character")
  
)



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);

#opt

#if (is.null(opt$file)){
#  print_help(opt_parser)
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}

#### TEST ANALYSIS
test<-0
if(test==1){
  
  #opt$prefix<-"TEST"
  #opt$directory<-"TWHO2.1/results"
  opt$metadata<-"TWHO_METADATA.txt"
  #opt$covariates<-"Study,Batch,Site,Collection_date,Type,Outbreak_Case"
  #opt$out<-"TWHO2.1/plots"
  #opt$tree<-"tree.nwk"
  #opt$covariates<-"Participant_ID,Site,Type,Outbreak_Case,Collection_date"
  #opt$textcovariates<-"Participant_ID,Collection_date"
  #opt$sort<-"Site,Type,Participant_ID"
  
  opt$prefix<-"TEST"
  opt$directory<-"TWHO_ALL/results_TG"
  opt$out<-"TWHO_ALL/results_TG"
  opt$tree<-"tree.nwk"
  opt$covariates<-"Participant_ID,Type,Outbreak_Case,Collection_date"
  opt$textcovariates<-"Participant_ID,Collection_date"
}


#opt$covariates<-"covariates_lineages.txt"
#opt$prefix<-"TWHO.lineage"

require(ggtree)
require(ggplot2)
suppressWarnings(require(patchwork))
suppressWarnings(require(reshape2))



# fix the alignment of the tree panel
#scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
#  scale_y_continuous(expand=expand, ...)
#}

scale_y_tree <- function(expand=expansion(0, 0.6), ...){
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

plot_variantsX<-function(a)
{
  a$pos = factor(alleles$pos)
  bases = c("A", "C", "G", "T", "N")
  a$alt_allele = factor(a$alt_allele, bases)
  ambiguous_bases = !(a$alt_allele %in% bases)
  a$alt_allele[ambiguous_bases] = "X"

  cols <- c("blue", "red", "green3", "purple3", "lightgrey")
  g <- ggplot(a, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    ylim(tip.order) +
    theme_bw() + 
    themes$variants + 
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  
  return(g) 
}

plot_variants<-function(a,s)
{
  a$pos = factor(a$pos)
  bases = c("A", "C", "G", "T", "N", "X")
  a$alt_allele = factor(a$alt_allele, bases)
  ambiguous_bases = !(a$alt_allele %in% bases)
  a$alt_allele[ambiguous_bases] = "X"
  
  cols <- c("blue", "red", "green3", "purple3", "lightgrey","black")
  g <- ggplot(a, aes(x=pos, y=name)) + 
    geom_tile(aes(fill=alt_allele), color="white") +
    ylim(s) +
    theme_bw() + 
    themes$variants + 
    scale_fill_manual(name="Variant", values=cols, drop=FALSE)
  g<-g + labs(title=opt$title)
  return(g) 
}



plot_covariate<-function(c,i)
{
  g<- ggplot(c,aes(x=1,y=ID)) + 
    geom_tile(aes_string(fill=i)) + 
    xlab(i) +
    themes$covariate +
    scale_fill_discrete(na.translate=FALSE)
  return(g)
}

plot_covariate_text<-function(c,i)
{
  g<- ggplot(c,aes(x=1,y=ID)) + 
    #geom_tile(aes(fill="white",colour="white")) + 
    geom_tile(colour="white",fill="white") + 
    geom_text(aes_string(label=i),size=1.6) +
    xlab(i) +
    themes$covariate +
    theme(legend.position = "none") +
    scale_fill_discrete(na.translate=FALSE)

  return(g)
}

get_tiporder<-function(tree){
  # get order of tips in the tree
  # from: https://groups.google.com/forum/#!topic/bioc-ggtree/LqRDK78m3U4
  d = fortify(tree)
  d = subset(d, isTip)
  tip.order = with(d, label[order(y, decreasing=F)])
  tip.order
}


get_covariates<-function(covariates_fn,sample.names,covariate.list){
  ### this will load the covariates file, restricting to only those that are requeste
  ### ID field corresponds to the sample name, in tip.order and allele table
  covariates<-NULL
  if(! is.null(covariates_fn)){
    cat("loading covariate data\n")
    #refdf<-data.frame(name="MN908947.3",location=NA,site=NA,batch=NA)
    covariates.all<-read.table(covariates_fn, header=T,as.is=T,fill=T)
    ### restrict to names
    covariates<-covariates.all[covariates.all$ID %in% sample.names,]
    
    ### add covariates that are missing from tip.order
    missing<-sample.names[! sample.names %in% covariates$ID]
    missing.df<-data.frame(matrix(NA,ncol=ncol(covariates),nrow=length(missing)))
    colnames(missing.df)<-colnames(covariates)
    missing.df$ID<-missing
    
    covariates<-rbind(covariates,missing.df)
    
    #rownames(covariates)<-covariates$ID
    #covariates<-covariates[,colnames(covariates) != "ID"]

    ### now restrict to requested covariates, if provided
    if(! is.null(covariate.list)){
      covariate.names<-c("ID",unlist(strsplit(covariate.list,",")))
      covariates<-covariates[,colnames(covariates) %in% covariate.names]

      missing<-covariate.names[!covariate.names %in% colnames(covariates)]

      if(length(missing)>0){
        msg<-paste("requested covariates not found in metadata:",paste(missing,collapse=","))
        warning(msg)
      }
      
      ### add the missing covariates, all as NA
      covariates[,missing]<-NA
      ### order the covariate columns
      covariates<-covariates[,covariate.names]
    }
  } 
  covariates
}


make_table<-function(a,c,o){
  m<-dcast(a,name~pos,value.var="alt_allele")
  m$name<-as.character(m$name)  ### remove factor
  m[is.na(m)]<-""
  
  ref<-unique(a[,c("pos","ref_allele")])
  ref<-ref[order(ref$pos),]
  refrow<-c("MN908947.3",as.character(ref$ref_allele))
  m<-rbind(refrow,m)
  rownames(m)<-m$name
  m<-m[rev(o),]
  
  rownames(c)<-c$name
  c<-c[rev(o),]
  c<-c[,-1]
  c[is.na(c)]<-""
  mm<-cbind(c,m)
  mm
}



#### MAIN
data_dir<-opt$directory
covariates_fn<-opt$metadata
alleles_fn<-opt$variants
prefix<-opt$prefix
plot_dir<-opt$out

if( (! is.null(opt$tree)) & (! is.null(opt$sort))){
  stop("detecting both a tree and a sort order.  Only one or the other should be indicated ")
}


### load the alleles, this forms the heatmap
cat("loading alleles\n")
alleles_path = paste(data_dir, alleles_fn, sep="/")
alleles <- read.table(alleles_path, header=T)
### get the sample names, 
### this must match a tree, if it is selected
### covariate data is limited ot these sample names
sample.names<-as.character(unique(alleles$name))
sample.names<-c(sample.names,opt$reference_name)
### get the covariates if provided, otherwise covariates will be null

covariates<-get_covariates(covariates_fn,sample.names,opt$covariates)


if(! is.null(opt$sort)){
  cat('applying sort order\n')
  sort.order<-unlist(strsplit(opt$sort,","))
  
  ## remove reference
  refdf<-covariates[covariates$ID==opt$reference_name,]
  covdf<-covariates[covariates$ID!=opt$reference_name,]
  covdf<-covdf[do.call("order",covdf[sort.order]),]
  
  covariates<-rbind(covdf,refdf)
  
  ### should remove the reference first, then add to the end
}
### reverse to put reference sequence at the bottom
sample.order<-rev(covariates$ID)



### get the tree plot, if this is provided
### sample order changes to order in the tree
if(! is.null(opt$tree)){
  cat("loading newick tree\n")
  tree_fn<-opt$tree 
  tree_path<-paste(data_dir, tree_fn, sep="/")
  tree <- read.tree(tree_path)
  
  ### get ordering from the tree
  sample.order<-get_tiporder(tree)
}


###
#### APPLY factor levels for ordering
alleles$name = factor(alleles$name, levels=sample.order)
covariates$ID=factor(covariates$ID,levels=sample.order)

### a list of plot panels, ordered from left to right
panels<-list()

### get the tree plo
if(! is.null(opt$tree)){
  panels$tree<-plot_tree(tree)
}

#covariate.ids<-colnames(covariates)
#covariate.ids<-covariate.ids[covariate.ids!="ID"]


#textcovariates<-unlist(strsplit(opt$textcovariates,","))

#for (id in covariate.ids){
#  if(id %in% textcovariates){
#      cat(id,"in text covariates\n")
#      panels[[id]]<-plot_covariate_text(covariates,id)
#  }else{
#      panels[[id]]<-plot_covariate(covariates,id)
#  }
#}


#### get the allele plot, a heatmap
#panels$variants<-plot_variants(alleles)   
panels$variants<-plot_variants(alleles,sample.order) 


# count number of samples, for scaling the plot
num_samples = length(sample.order)
cat("number of samples for this heatmap is ",num_samples,"\n")
num_sites = length(unique(alleles$pos))
cat("number of genomic positions for this heatmap is ",num_sites,"\n")
#num_covariates<-length(covariate.ids)
#cat("number of covariates for this heatmap is ",num_covariates,"\n")



### set the panel width, 1 for the tree, 0.1 for each covariate and 4 for the variants
### the proportion for the variant panel should adjust based on the number of variants
#panel.widths<-c(1,rep(0.1,length(covariate.ids)),4)
#covariate_panel_widths<-rep(0.1,num_covariates)
#covariate_panel_widths[covariate.ids %in% textcovariates]<-0.3

panel.widths<-NULL
image.width<-0
if(! is.null(opt$tree)){
  panel.widths<-c(panel.widths,1)
  image.width<-1
}
#panel.widths<-c(panel.widths,covariate_panel_widths)
panel.widths<-c(panel.widths,0.1*num_sites)

#image.width<-image.width + sum(covariate_panel_widths)
image.width<-image.width + (0.05*num_sites)*4

image.height<-0.25 * num_samples
image.height<-max(image.height,5)


if(! is.null(opt$height)){
  image.height<-as.numeric(opt$height)
}
if(! is.null(opt$width)){
  image.width=as.numeric(opt$width)
}


cat("calculated image width ",image.width,"\n")
cat("calculated image height ",image.height,"\n")


p<-wrap_plots(panels) + plot_layout(widths = panel.widths, guides="collect")


export.table<-make_table(alleles,covariates,sample.order)
### save data table
fn<-paste(plot_dir,"/",prefix,"_variant_map.tsv",sep="")
write.table(export.table,file=fn,quote=F,sep="\t",row.names=F,col.names=T)
### save RData
fn<-paste(plot_dir,"/",prefix,"_variant_map.Rdata",sep="")
save.image(file=fn)
## save the plot    
fn = paste(plot_dir,"/",prefix,"_variant_map.pdf",sep="")
## wdith shoud be a factor of the number of variants
ggsave(filename=fn,plot=p,height=image.height, width=image.width, limitsize=F)
