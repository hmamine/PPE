#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "vegan", "metagenomeSeq","PMCMRplus")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c("#155832","#ffa500","#00b2b2","#4f7b95")
shape.1=c(17,19)

theme_new <- theme (
	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	)
	
theme_new1 <- theme (
	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.text.x=element_blank(),
	)
	
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 
meth1=("PCoA")
meth2=("CAP")
dist1=("bray")

set.seed(123)

# upload and prepare phyloseq objects***
	mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	TREE=read.tree("tree.nwk")

	physeq=phyloseq(OTU, TAXA, SD) 
	physeq1=subset_samples( physeq, Treat %in% c("AbxEa","control","Ea","PPE") )
	physeq.B=subset_taxa( physeq1, !Genus %in% c("Mitochondria","Chloroplast") )
	Taxa=tax_table(physeq.B)
	Sd=sample_data(physeq.B)
#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.B ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), Taxa, Sd ,TREE )

#computing unconstrained PCoA using Bray-Curtis distances

	dist_BC=distance( physeq_norm, dist1 )
	pcoa_BC=ordinate( physeq_norm, meth1 , dist_BC )
	pBC<-plot_ordination( physeq_norm, pcoa_BC, color="Treat", shape="Type")	
	pBC$layers<-pBC$layers[-1]
	
	p0=pBC+geom_point(size=4)+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)
	legend <- g_legend(p0) 
	
	p1=pBC+geom_point(size=4)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)

#	pdf("pcoa_plot.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(p3, nrow=2, ncol=2)
#	dev.off()
	
	
	
	
	
	
	
