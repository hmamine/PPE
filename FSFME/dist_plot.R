#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime", "reshape2")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

list.color<-c("#b02b76","#ffa500","#00b2b2")

DT<-data.table(rn=c("S1","S2","S3","S4","S5"),  key="rn")
q<-list()

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_blank()
	)
dist=( "bray" )

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

	physeq_B=subset_taxa( physeq, Family  != "mitochondria" )

	physeq.subset=subset_samples( physeq_B, treatment %in% c("M","I","A","L") )

LIST1<-list ("ChSp","Obx")
LIST2<-list ("T4","T6","T8")
BC<-list()
for (j in LIST1)
{
print( j )
SUB.j=c( j )
physeq_j=subset_samples( physeq.subset, cultivar %in% SUB.j )
	
	for( i in LIST2 ) 
	{
	print( i )
	SUB.i=c( i )
	physeq_i=subset_samples( physeq_j, time %in% SUB.i )

	#trimming OTUs with low occurance
	otu.I=otu_table(physeq_i)
	taxa.I=tax_table(physeq_i)
	sd.I=sample_data(physeq_i)

	flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__","sc__","ot_")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
TAXA.filt[,"Phylum"] <- ifelse(! TAXA.filt[,"Phylum"] %in% ColPhylum, "Other", ifelse(TAXA.filt[,"Phylum"] == "Bacteroidetes", "Bacteroidetes", ifelse(TAXA.filt[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (TAXA.filt[,"Phylum"] == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter.I=phyloseq(NewOTU, NewTaxa, sd.I)

	#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.filter.I ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), sd.I, taxa.I, TREE )
	
	#computing distances
	DT.sd=data.table( as( sample_data( physeq_norm ), "data.frame" ),keep.rownames=T, key="rn" )
	dist_BC=distance( physeq_norm, dist )		
	
	M.dist=melt(as.matrix( dist_BC ))

	M.dist = M.dist %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)

	DT.dist=data.table(M.dist,  key=c("Var1","Var2"))
	DT.merge=merge( DT.dist, DT.sd, by.x="Var2", by.y="rn" )
	DT.merge=merge( DT.merge, DT.sd, by.x="Var1", by.y="rn" )

	DT1=DT.merge[treatment2.x == "M" & treatment2.y != "M"]
	
	BC[[j]][[i]]$I<-DT1[DT1$treatment.y == "I"]$value
	BC[[j]][[i]]$A<-DT1[DT1$treatment.y == "A"]$value
	BC[[j]][[i]]$L<-DT1[DT1$treatment.y == "L"]$value
	
}
}

	D.Obx<-rbind(BC$Obx$T4$I, BC$Obx$T4$A,BC$Obx$T4$L, BC$Obx$T6$I, BC$Obx$T6$A,BC$Obx$T6$L, BC$Obx$T8$I, BC$Obx$T8$A, BC$Obx$T8$L)
	rownames(D.Obx)<-paste("S", rep (1:9), sep ="")
	colnames(D.Obx)<-paste("D", rep (1:81), sep ="")
	D.Obx=as.data.frame(D.Obx)
	D.Obx$time=c("T4","T4","T4","T6","T6","T6","T8","T8","T8")
	D.Obx$treat=c("I","A","L","I","A","L","I","A","L")

	DT_Obx<-data.table(melt(D.Obx))
	DT_Obx$Concat<-paste(DT_Obx$time,DT_Obx$treat, sep="")

	pp2=ggplot(DT_Obx, aes(x=Concat, y=value))
	pObx=pp2+geom_violin(aes(fill=treat), color="black")+
	geom_smooth(aes(group=treat, linetype=treat),colour="Black",method=glm,se=F)+
	geom_boxplot(width=0.1, outlier.shape = 1, outlier.size= 1)+
	theme_bw()+ theme_new +
	scale_y_continuous(limits=c(.3,1), breaks=c(.3,.5,.7,.9))+
	scale_fill_manual(values=list.color)+
	scale_colour_manual(values=list.color)+
	ylab("BC distances to mock")

	DT_ObxT4<-DT_Obx[DT_Obx$time == "T4",]
	Obx.T4<-aov(data=DT_ObxT4, value ~ treat)
	print(summary(Obx.T4))
	print(HSD.test(Obx.T4,trt="treat"))

	DT_ObxT6<-DT_Obx[DT_Obx$time == "T6",]
	Obx.T6<-aov(data=DT_ObxT6, value ~ treat)
	print(summary(Obx.T6))
	print(HSD.test(Obx.T6,trt="treat"))

	DT_ObxT8<-DT_Obx[DT_Obx$time == "T8",]
	Obx.T8<-aov(data=DT_ObxT8, value ~ treat)
	print(summary(Obx.T8))
	print(HSD.test(Obx.T8,trt="treat"))


D.ChSp<-rbind(BC$ChSp$T4$I, BC$ChSp$T4$A, BC$ChSp$T4$L, BC$ChSp$T6$I, BC$ChSp$T6$A,BC$ChSp$T6$L, BC$ChSp$T8$I, BC$ChSp$T8$A, BC$ChSp$T8$L)
	rownames(D.ChSp)<-paste("S", rep (1:9), sep ="")
	colnames(D.ChSp)<-paste("D", rep (1:81), sep ="")
	D.ChSp=as.data.frame(D.ChSp)
	D.ChSp$time=c("T4","T4","T4","T6","T6","T6","T8","T8","T8")
	D.ChSp$treat=c("I","A","L","I","A","L","I","A","L")

	DT_ChSp<-data.table(melt(D.ChSp))
	DT_ChSp$Concat<-paste(DT_ChSp$time,DT_ChSp$treat, sep="")

	pp1=ggplot(DT_ChSp, aes(x=Concat, y=value))
	pChSp=pp1+geom_violin(aes(fill=treat), color="black")+
	geom_smooth(aes(group=treat, linetype=treat),colour="Black",method=glm,se=F)+
	scale_y_continuous(limits=c(.3,1), breaks=c(.3,.5,.7,.9))+
	geom_boxplot(width=0.1, outlier.shape = 1, outlier.size= 1)+	
	theme_bw()+ theme_new +
	scale_fill_manual(values=list.color)+
	ylab("BC distances to mock")

	DT_ChSpT4<-DT_ChSp[DT_ChSp$time == "T4",]
	ChSp.T4<-aov(data=DT_ChSpT4, value ~ treat)
	print(summary(ChSp.T4))

	DT_ChSpT6<-DT_ChSp[DT_ChSp$time == "T6",]
	ChSp.T6<-aov(data=DT_ChSpT6, value ~ treat)
	print(summary(ChSp.T6))
	print(HSD.test(ChSp.T6,trt="treat"))

	DT_ChSpT8<-DT_ChSp[DT_ChSp$time == "T8",]
	ChSp.T8<-aov(data=DT_ChSpT8, value ~ treat)
	print(summary(ChSp.T8))


#	pdf("BC_distances_plots.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(pObx, pChSp, nrow=2, ncol=2)
#	dev.off()










	

























