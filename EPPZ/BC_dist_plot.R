#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table","agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime","picante","PMCMRplus")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
color_palette<-c("#ffa500","#00b2b2","#4f7b95")



theme_new <- theme (
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.ticks.x=element_blank(),
	axis.text.x = element_blank(),
#	axis.text.y = element_blank(),
	legend.position="none",
#	axis.text.x = element_text(angle=90, vjust=1, size=9, color="black"),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
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

	physeq=phyloseq(OTU, TAXA, SD) 
	physeq1=subset_samples( physeq, Treat %in% c("control","Ea","PPE") )
	physeq.B=subset_taxa( physeq1, !Genus %in% c("Mitochondria","Chloroplast") )
	Taxa=tax_table(physeq.B)
	Sd=sample_data(physeq.B)
#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.B ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), Taxa, Sd ,TREE )

#computing BC distances for 16S rRNA
	physeq.R=subset_samples( physeq_norm, Type == "RNA")
	DT.sd.R=data.table(as(sample_data(physeq.R), "data.frame"),keep.rownames=T, key="rn")
	dist.BC.R=distance(physeq.R, dist )
	melt.dist.BC.R=reshape2::melt(as.matrix( dist.BC.R ))

	DT.dist.R = melt.dist.BC.R %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)
 	
   	DT.dist.R=data.table(DT.dist.R ,  key=c("Var1","Var2"))
   	DT.R=merge( DT.dist.R, DT.sd.R, by.x="Var2", by.y="rn" )
	DT.R=merge( DT.R, DT.sd.R, by.x="Var1", by.y="rn" )	
	
	DT.R=DT.R[DT.R$Treat.x=="control" & DT.R$Treat.y!="control"]
	pp1=ggplot(DT.R, aes(x=Concat.y, y=value))

	p1=pp1+geom_violin(aes(fill=Treat.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette)+
	theme_bw()+ theme_new + ylab("BC distances - RNA")
	print(p1)
	
	
#computing BC distances for 16S rDNA
	physeq.D=subset_samples( physeq_norm, Type == "DNA")
	DT.sd.D=data.table(as(sample_data(physeq.D), "data.frame"),keep.rownames=T, key="rn")
	dist.BC.D=distance(physeq.D, dist )
	melt.dist.BC.D=reshape2::melt(as.matrix( dist.BC.D ))

	DT.dist.D = melt.dist.BC.D %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)
 	
   	DT.dist.D=data.table(DT.dist.D ,  key=c("Var1","Var2"))
   	DT.D=merge( DT.dist.D, DT.sd.D, by.x="Var2", by.y="rn" )
	DT.D=merge( DT.D, DT.sd.D, by.x="Var1", by.y="rn" )	
	
	DT.D=DT.D[DT.D$Treat.x=="control" & DT.D$Treat.y!="control"]
	pp2=ggplot(DT.D, aes(x=Concat.y, y=value))

	p2=pp2+geom_violin(aes(fill=Treat.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette)+
	theme_bw()+ theme_new + ylab("BC distances")
	print(p2)


	
	DT=rbind(DT.R,DT.D)
	pp3=ggplot(DT, aes(x=Concat.y, y=value))
	Ord=c("DNA.Ea","RNA.Ea","DNA.PPE","RNA.PPE")
	pp3$data$Concat.y <- ordered(pp3$data$Concat.y, levels=Ord )
	p3=pp3+geom_violin(aes(fill=Treat.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette)+
	theme_bw()+ theme_new + ylab("BC distances")
	print(p2)
	
	
#	pdf("BC_distances.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(p3, nrow=2, ncol=2)
#	dev.off()
	
	kruskal.test(data=p3$data, value ~ Concat.y)
	kwAllPairsConoverTest(value ~ as.factor(Concat.y) , data=p3$data,  p.adjust.method="BH" )
#	aov.out<-aov(data=p3$data, value ~ Concat.y)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Concat.y"))






















	

