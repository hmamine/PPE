#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table","agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c("#cc2a36","#166590","#11902d")

theme_new <- theme (
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.ticks.x=element_blank(),
	axis.text.x = element_blank(),
#	axis.text.y = element_blank(),
	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
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

#computing distance to centroid
	DT.sd=data.table( as( sample_data( physeq_norm ), "data.frame" ),keep.rownames=T, key="rn" )
	dist_i=distance( physeq_norm, dist )		
	disp <- with( DT.sd, betadisper( dist_i, concat2, type="centroid" ) )	
#	print(disp)
#	print(anova(disp))
#	print(TukeyHSD(disp))

	DT.data=data.table(as.data.frame(disp$distances), keep.rownames=T, key="rn")
	setnames(DT.data, "disp$distances","value")
	DT.sd=data.table( as( Sd, "data.frame" ),keep.rownames=T, key="rn" )
	DT=DT.data[DT.sd]

	
	pp1=ggplot(DT, aes(x=concat2, y=value))
	

	neworder=c("untreat","2D2dpi_H2O","2D2dpi_Ea","2D4dpi_H2O","2D4dpi_Ea","4D2dpi_H2O","4D2dpi_Ea")
	pp1$data$concat2 <- ordered(pp1$data$concat2, levels=neworder )
	
	p1=pp1+geom_boxplot(aes(fill=Treat),width=0.75)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_fill_manual(values=color.1)+
	ylab("BC distances to centroid")
	
#	print(kruskal.test(data=p1$data[p1$data$concat1=="2D2dpi",], value ~ Treat))
#	DT1=p1$data[p1$data$concat1 == "2D2dpi",]
#	summary(out.B<-aov(data=p1$data[p1$data$concat1=="2D2dpi",], value ~ as.factor(Treat)))
#	wilcox.test( DT1$value[DT1$Treat == "Ea"], DT1$value[DT1$Treat == "H20"] , correct=T, exact=T)
#	print(kruskal.test(data=p1$data[p1$data$concat1=="2D4dpi",], value ~ Treat))
#	print(kruskal.test(data=p1$data[p1$data$concat1=="4D2dpi",], value ~ Treat))
	
	




	dist_BC=distance( physeq_norm, dist )		
	
	M.dist=melt(as.matrix( dist_BC ))
	
	M.dist = M.dist %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)
   	
   	DT.dist=data.table(M.dist,  key=c("Var1","Var2"))
   	DT.merge=merge( DT.dist, DT.sd, by.x="Var2", by.y="rn" )
	DT.merge=merge( DT.merge, DT.sd, by.x="Var1", by.y="rn" )

	DT1=DT.merge[concat1.x == "2D2dpi" & concat1.y =="2D2dpi"]
	DT1=DT1[Treat.x == "H20"]
#	DT1=DT1[,c(7,3)]

	DT2=DT.merge[concat1.x == "2D4dpi" & concat1.y =="2D4dpi"]
	DT2=DT2[Treat.x == "H20"]
#	DT2=DT2[,c(7,3)]
	
	DT3=DT.merge[concat1.x == "4D2dpi" & concat1.y =="4D2dpi"]
	DT3=DT3[Treat.x == "H20" ]
#	DT3=DT3[,c(7,3)]
	
	DT=rbind(DT1,DT2,DT3)
	 
	p2=ggplot(DT, aes(x=concat1.x, y=value, fill=Treat.y))+
	geom_violin()+geom_boxplot(colorr="white",width=0.1,outlier.size =0.5)+
	theme_bw()+ theme_new + ylab("BC distances") +
	scale_fill_manual(values=color.1)

	
#	aov.out<-aov(data=DT, value ~ concat2.x)
#	print(summary(aov.out))
#	print(HSD.test(aov.out,trt="concat2.x"))
#	kruskal.test(data=DT, value ~ concat2.x)



#	pdf("BC_distances_plot.pdf",paper=,useDingbats=FALSE)
	gridExtra::grid.arrange(p1,p2, nrow=2, ncol=2)
#	dev.off()

#Permanova test for each dpi to mock inoculated samples	
	physeq1=subset_samples( physeq_norm, concat1 == "2D2dpi")
	dist1=distance(physeq1, dist)
	sd1=data.table( as( sample_data(physeq1), "data.frame" ),keep.rownames=T, key="rn" )
	ad.out.1<-with( sd1, adonis ( dist1 ~ Treat ) )
	print(ad.out.1)
	
	physeq2=subset_samples( physeq_norm, concat1 == "2D4dpi")
	dist2=distance(physeq2, dist)
	sd2=data.table( as( sample_data(physeq2), "data.frame" ),keep.rownames=T, key="rn" )
	ad.out.2<-with( sd2, adonis ( dist2 ~ Treat ) )
	print(ad.out.2)
	
	physeq3=subset_samples( physeq_norm, concat1 == "4D2dpi")
	dist3=distance(physeq3, dist)
	sd3=data.table( as( sample_data(physeq3), "data.frame" ),keep.rownames=T, key="rn" )
	ad.out.3<-with( sd3, adonis ( dist3 ~ Treat ) )
	print(ad.out.3)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	



