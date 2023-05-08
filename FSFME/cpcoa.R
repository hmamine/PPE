#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape","metagenomeSeq","PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c("#cc2a36","#166590","#11902d")
shape.1=c(17,19,15,13)

theme_new <- theme (
#	axis.title.y=element_blank(),
#	axis.title.x=element_blank(),
#	axis.text.x = element_blank(),
#	axis.text.y = element_blank(),
#	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	)
meth1=("PCoA")
meth2=("CAP")
dist1=("bray")

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

	physeq= phyloseq(OTU, TAXA, SD,TREE)

	physeq.B=subset_samples( physeq, Treat != "Negative" )
	physeq.BM=subset_taxa( physeq.B, !Genus %in% c("Mitochondria","Chloroplast") )
	physeq.BT=subset_taxa( physeq.BM, Select != "Prune" )
	Taxa=tax_table(physeq.BT)
	Sd=sample_data(physeq.BT)

#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.BT ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), Taxa, Sd ,TREE )

#computing Bray-Curtis distances
	dist_BC=distance( physeq_norm, dist1 )

#computing constrained PCoA	
	cpcoa_BC=ordinate( physeq_norm, meth2 , dist_BC, ~concat2 )
	cpBC<-plot_ordination( physeq_norm, cpcoa_BC, color="Treat", shape="concat1")
	cpBC$layers<-cpBC$layers[-1]
	cpBC$data <- merge(cpBC$data,aggregate(cbind(mean.CAP1=CAP1,mean.CAP2=CAP2)~concat2,cpBC$data,mean),by="concat2")
	
	p1=cpBC+geom_point(size=3)+theme_bw()+theme_new+
	geom_segment(aes(x=mean.CAP1, y=mean.CAP2, xend=CAP1, yend=CAP2),size=0.3)+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)
	
#print(anova.cca(cpcoa_BC.L, by="term", perm.max=1000))
#DT.L=data.table( as( sample_data( physeq_norm.L ), "data.frame" ),keep.rownames=T, key="rn" )
#with( DT.L, adonis ( dist_BC.L ~ Geno ) )



#	pdf("pcoa_plot.pdf",useDingbats=FALSE)
	print(p1)
#	dev.off()

