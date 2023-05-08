#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "ape", "metagenomeSeq","PMCMR","UpSetR", "venn","seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")

Shape.1=c(21,19)

theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	legend.position="none"
	)

theme_new2 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	legend.position="none"
	)

Manhat<-function(mapping){
		p<- ggplot( DT.taxa , mapping )
		p+geom_jitter( size=2 )+theme_bw( )+
		geom_hline(yintercept = -log2(0.05), linetype="longdash")+
		scale_color_manual(values=Palette.phylum)+
		scale_shape_manual(values=Shape.1)+
		coord_flip()+theme_new1
	}

Volcan<-function(mapping){
		p<- ggplot( DT.taxa , mapping )
		p+geom_jitter( size=2 )+ theme_bw( )+
		scale_color_manual(values=Palette.phylum)+
		geom_hline(yintercept = -log2(0.05), linetype="longdash")+
		geom_vline(xintercept = c(-1,0,1), linetype="longdash")+
		theme_new1
	}

res_ad1<-vector( "list" )
res_ad2<-vector( "list" )
res_ad3<-vector( "list" )

#log2 fold change
	folch=1

#Adjusted p-value
	alpha=0.05

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




physeq.2D=subset_samples( physeq.BT, concat1 != "4D2dpi" )

#trimming OTUs with low occurance
otu.I=otu_table(physeq.2D)
taxa.I=tax_table(physeq.2D)
sd.I=sample_data(physeq.2D)
flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
otu.I.filtered=flt$mat
taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
dummyTaxonomy=c("D__dummy","Ph__","Cl__","Or__","Fa__","Ge__","Sp__","Se__")
TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
#TAXA.filt[,"Phylum"] <- ifelse(! TAXA.filt[,"Phylum"] %in% ColPhylum, "Other", ifelse(TAXA.filt[,"Phylum"] == "Bacteroidetes", "Bacteroidetes", ifelse(TAXA.filt[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (TAXA.filt[,"Phylum"] == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
NewTaxa=tax_table(TAXA.filt)

physeq.flt=phyloseq(NewOTU, NewTaxa, sd.I)

physeq.flt=tax_glom(physeq.flt, "Family") 

otu.i=otu_table(physeq.flt)
taxa.i=tax_table(physeq.flt)
DT.taxa=data.table(taxa.i, keep.rownames=T, key="rn" )
sd.i=sample_data(physeq.flt)

Time<-list ( "2D2dpi","2D4dpi")
for( i in Time ) 
{
	print( i )
#physeq.I=subset_samples( physeq.subset.Loc, Time %in% subset.I )
	physeq.I=subset_samples( physeq.flt, concat1 == i )
	

#computing differentially abundant OTUs 
	m=as( otu_table( physeq.I ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.I ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.I ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Treat
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Treat, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( TreatH20	-	TreatEa, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")
	colnames(DT.zig) <- paste(colnames(DT.zig), i , sep = ".")

#		DT_1=data.table( DT_1, keep.rownames=T, key="rn" )	
#		DT_1p=data.table(DT_1p, keep.rownames=T, key="rn" )
#		setnames( DT_1, "Treatmentzt - Treatmentmock", paste0( "fc.ZtvsM.",i ) )	
#		setnames( DT_1p, "Treatmentzt - Treatmentmock", paste0( "p.ZtvsM.",i ) )
#		DT_1=merge( DT_1, DT_1p , by="rn" )

	DT.taxa=DT.taxa[ DT.zig, ]
}
break()
	ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")

	DT.taxa$ColPhylum <- ifelse(!DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Bacteroidetes", "Bacteroidetes", ifelse(DT.taxa$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.taxa$Phylum == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))

	DT.taxa$Sh.T4 <- ifelse ( DT.taxa$logFC.T4 > 0, "Enriched", "Depleted")
	DT.taxa$Sh.T6 <- ifelse ( DT.taxa$logFC.T6 > 0, "Enriched", "Depleted")
	DT.taxa$Sh.T8 <- ifelse ( DT.taxa$logFC.T8 > 0, "Enriched", "Depleted")
	
	res_ad1[[j]]<-DT.taxa
	
	p1<-Manhat( aes( y=-log2(adj.P.Val.T4), x=ColPhylum, color = ColPhylum, shape=Sh.T4 ))
	p2<-Manhat( aes( y=-log2(adj.P.Val.T6), x=ColPhylum, color = ColPhylum, shape=Sh.T6 ))
	p3<-Manhat( aes( y=-log2(adj.P.Val.T8), x=ColPhylum, color = ColPhylum, shape=Sh.T8 ))
#	gridExtra::grid.arrange(p1, p2, p3, ncol=3)
	
	p4<-Volcan( aes( y=-log2(adj.P.Val.T4), x=logFC.T4, color = ColPhylum ))
	p5<-Volcan( aes( y=-log2(adj.P.Val.T6), x=logFC.T6, color = ColPhylum ))
	p6<-Volcan( aes( y=-log2(adj.P.Val.T8), x=logFC.T8, color = ColPhylum ))
#	gridExtra::grid.arrange(p4, p5, p6, ncol=3)
	
	T4<-DT.taxa[ which(DT.taxa$adj.P.Val.T4 < alpha & abs( DT.taxa$logFC.T4 ) > folch ), ]$rn
	T6<-DT.taxa[ which(DT.taxa$adj.P.Val.T6 < alpha & abs( DT.taxa$logFC.T6 ) > folch ), ]$rn
	T8<-DT.taxa[ which(DT.taxa$adj.P.Val.T8 < alpha & abs( DT.taxa$logFC.T8 ) > folch ), ]$rn
	
#	pdf(paste0("venn_diagram_",j,".pdf"), useDingbats=FALSE)
	VennPlot<- venn(list(T4, T6, T8), zcolor="style")
#	dev.off()
	InterSec<-attr(VennPlot, "intersection")
	DT.InterSec<- data.table(DT.taxa[DT.taxa$rn %in% InterSec$'A:B:C',], key="rn")
	res_ad2[[j]]<-DT.InterSec

	Uniq<-unique(c(T4,T6,T8), fromLast = TRUE)
	DT.Uniq<-DT.taxa[ which(DT.taxa$rn %in% Uniq ), ]
	res_ad3[[j]]<-DT.Uniq



DT.I<-res_ad1[[1]]
T4.I<-DT.taxa[ which(DT.I$adj.P.Val.T4 < alpha & abs( DT.I$logFC.T4 ) > folch ), ]$rn
T6.I<-DT.taxa[ which(DT.I$adj.P.Val.T6 < alpha & abs( DT.I$logFC.T6 ) > folch ), ]$rn
T8.I<-DT.taxa[ which(DT.I$adj.P.Val.T8 < alpha & abs( DT.I$logFC.T8 ) > folch ), ]$rn


DT.A<-res_ad1[[2]]
T4.A<-DT.taxa[ which(DT.A$adj.P.Val.T4 < alpha & abs( DT.A$logFC.T4 ) > folch ), ]$rn
T6.A<-DT.taxa[ which(DT.A$adj.P.Val.T6 < alpha & abs( DT.A$logFC.T6 ) > folch ), ]$rn
T8.A<-DT.taxa[ which(DT.A$adj.P.Val.T8 < alpha & abs( DT.A$logFC.T8 ) > folch ), ]$rn

#pdf(paste0("venn_diagram_4dpi_",l,".pdf"), useDingbats=FALSE)
VennPlotT4<-venn(list(T4.I, T4.A), zcolor="style")
dev.off()
InterSecT4<-attr(VennPlotT4, "intersection")
DT.InterSecT4<- data.table(DT.taxa[DT.taxa$rn %in% InterSecT4$'A:B',], key="rn")

#pdf(paste0("venn_diagram_6dpi_",l,".pdf"), useDingbats=FALSE)
VennPlotT6<-venn(list(T6.I, T6.A), zcolor="style")
dev.off()
InterSecT6<-attr(VennPlotT6, "intersection")
DT.InterSecT6<- data.table(DT.taxa[DT.taxa$rn %in% InterSecT6$'A:B',], key="rn")

#pdf(paste0("venn_diagram_8dpi_",l,".pdf"), useDingbats=FALSE)
VennPlotT8<-venn(list(T8.I, T8.A), zcolor="style")
dev.off()
InterSecT8<-attr(VennPlotT8, "intersection")
DT.InterSecT8<- data.table(DT.taxa[DT.taxa$rn %in% InterSecT8$'A:B',], key="rn")

print("done")


T4.cv<-DT.InterSecT4$rn
T6.cv<-DT.InterSecT6$rn
T8.cv<-DT.InterSecT8$rn

venn(list(T4.cv, T6.cv, T8.cv), zcolor="style")





#T4<-IObx[ which(IObx$adj.P.Val.T4 < alpha & abs( IObxa$logFC.T4 ) > folch ), ]$rn
#T6<-IObx[ which(IObx$adj.P.Val.T6 < alpha & abs( IObx$logFC.T6 ) > folch ), ]$rn
#T8<-IObx[ which(DT.taxa$adj.P.Val.T8 < alpha & abs( DT.taxa$logFC.T8 ) > folch ), ]$rn














