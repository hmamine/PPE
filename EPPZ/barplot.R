#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "metagenomeSeq", "magrittr","RColorBrewer","PMCMRplus")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
pl<-list()
color.1<-c("#b02b76","#155832","#ffa500","#00b2b2","#4f7b95")
color.2<-c("#ffa500","#6c6b20","#b36200")
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position ="none", 
	axis.ticks.x =  element_blank(),
	panel.background = element_blank(),
	panel.grid =  element_blank(),
	panel.border = element_blank(),
#	axis.text.x = element_text(angle=90, vjust=1),
	axis.title.x=element_blank(),
	axis.text.x=element_blank(),
	axis.title.y=element_blank(),
	)
theme_new1 <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
#	legend.position ="none", 
	axis.ticks.x =  element_blank(),
	panel.background = element_blank(),
	panel.grid =  element_blank(),
	axis.title.x=element_blank(),
#	axis.text.x=element_blank(),
#	axis.title.y=element_blank(),
	)
	
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

###transform counts to rel. abund.	
	physeq.RA = transform_sample_counts(physeq.B, function(x) x/sum(x))

###agglomerate at the family level
#	physeq.agg=tax_glom(physeq.RA, "Genus") 

#Agglomeration by sample type
#	physeq_agg=merge_samples(physeq.RA, "concat1") 
#	physeq.Prot=subset_taxa( physeq_agg, Phylum == "Proteobacteria" )


	m=as( otu_table( physeq.RA ), "matrix" )
	DT_m <- m %>% 
	melt(id.vars = "variable", value.name="relabund")
	DT.m=data.table(DT_m, keep.rownames=F, key="Var2") 
	
	DT.taxa=data.table( as(TAXA,"matrix"), keep.rownames=T, key="rn")
	DT.SD=data.table( as(SD,"data.frame"), keep.rownames=T, key="rn")
	
	DT=merge(DT.m, DT.taxa, by.x="Var1",by.y="rn")
	DT=merge(DT, DT.SD, by.x="Var2",by.y="rn")

	ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
	DT$ColPhylum <- ifelse(!DT$Phylum %in% ColPhylum, "Other", ifelse(DT$Phylum == "Bacteroidota", "Bacteroidetes", ifelse(DT$Phylum == "Firmicutes", "Firmicutes", ifelse (DT$Phylum == "Actinobacteriota", "Actinobacteria", "Proteobacteria"))))
	
	DT$Col<-ifelse(DT$relabund > 0.01, paste0(DT$Genus,DT$Var1), "Others")

	getPalette = colorRampPalette(brewer.pal(12, "Paired"))
	Len=length(unique(DT$Col))+1
	Palette.Genus<-getPalette(Len)

	pp=ggplot(DT, aes(x=Concat, y=relabund, fill=Col, color=Col ))

	Ord1<-c("Rcontrol","RAbxEa","REa","RPPE","Dcontrol","DAbxEa","DEa","DPPE")  
	pp$data$Concat <- ordered(pp$data$Concat, levels=Ord1 )
	
	p1= pp + geom_bar(position="fill", stat="identity")+
	scale_fill_manual(values=Palette.Genus)+
	scale_color_manual(values=Palette.Genus) 
	legend <- g_legend(p1) 
	
	p1= pp + geom_bar(position="fill", stat="identity") + theme_bw() +
	theme_new +
	scale_fill_manual(values=Palette.Genus)+
	scale_color_manual(values=Palette.Genus) +
	ylab("Relative Abundance")

#	pdf("bar_plot.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(p1, legend, nrow=2, ncol=2)
#	print(p1)
#	dev.off()

	DT1<-DT[DT$relabund > 0.01]
	
	
	Pa=DT1[DT1$Var1 == "4434d904522c3f016aa715d906835ac9",]
	pp3=ggplot(Pa, aes(x=Concat, y=relabund, color=Treat))
	pp3$data$Concat <- ordered(pp3$data$Concat, levels=Ord1 )
	pPa=pp3+ geom_boxplot()+ geom_point(size=2)+
	theme_bw()+ theme_new1 +
	scale_color_manual(values=color.1) +
	ylab("Relative Abundance of Pantoea ASVs")
	
	kruskal.test(data=pPa$data[pPa$data$Type=="D"], relabund ~ Treat)
	kwAllPairsConoverTest(relabund ~ as.factor(Treat) , data=pPa$data[pPa$data$Type =="D"],  p.adjust.method="BH" )
	
	kruskal.test(data=pPa$data[pPa$data$Type=="R"], relabund ~ Treat)
	kwAllPairsConoverTest(relabund ~ as.factor(Treat) , data=pPa$data[pPa$data$Type =="R"],  p.adjust.method="BH" )
	
	
	Ps=DT1[DT1$Var1 == "e1560e5a84cf55cda6e9e1a1ad09eb07",]
	pp4=ggplot(Ps, aes(x=Concat, y=relabund, color=Treat))
	pp4$data$Concat <- ordered(pp4$data$Concat, levels=Ord1 )
	pPs=pp4+ geom_boxplot()+ geom_point(size=2)+
	theme_bw()+ theme_new1 +
	scale_color_manual(values=color.1) +
	ylab("Relative Abundance of Pseudomonas ASVs")
	
	kruskal.test(data=pPs$data[pPs$data$Type=="D"], relabund ~ Treat)
	
	kruskal.test(data=pPs$data[pPs$data$Type=="R"], relabund ~ Treat)
	kwAllPairsConoverTest(relabund ~ as.factor(Treat) , data=pPs$data[pPs$data$Type =="R"],  p.adjust.method="BH" )
	
	
	Ea=DT1[DT1$Var1 == "9dc299594ffefeeb80611678432c6d94",]
	pp5=ggplot(Ea, aes(x=Concat, y=relabund, color=Treat))
	pp5$data$Concat <- ordered(pp5$data$Concat, levels=Ord1 )
	pEa=pp5+ geom_boxplot()+ geom_point(size=2)+
	theme_bw()+ theme_new1 +
	scale_color_manual(values=color.1) +
	ylab("Relative Abundance of Erwinia ASVs")
	
	kruskal.test(data=pEa$data[pEa$data$Type=="D"], relabund ~ Treat)
	kwAllPairsConoverTest(relabund ~ as.factor(Treat) , data=pEa$data[pEa$data$Type =="D"],  p.adjust.method="BH" )
	
	kruskal.test(data=pEa$data[pEa$data$Type=="R"], relabund ~ Treat)
	kwAllPairsConoverTest(relabund ~ as.factor(Treat) , data=pEa$data[pEa$data$Type =="R"],  p.adjust.method="BH" )
	
#	pdf("box_plot_ASVs.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(pPa,pPs,pEa, nrow=2, ncol=2)
#	dev.off()
	
	
	Ea.m.DCtl=mean(Ea[Ea$Concat =="Dcontrol"]$relabund)
	Ea.m.DEa=mean(Ea[Ea$Concat =="DEa"]$relabund)
	Ea.m.DPPE=mean(Ea[Ea$Concat =="DPPE"]$relabund)
	
	Ea.m.RCtl=mean(Ea[Ea$Concat =="Rcontrol"]$relabund)
	Ea.m.REa=mean(Ea[Ea$Concat =="REa"]$relabund)
	Ea.m.RPPE=mean(Ea[Ea$Concat =="RPPE"]$relabund)
	
	print(paste0("Ea-Ctl= ",Ea.m.RCtl/Ea.m.DCtl))
	print(paste0("Ea-Ea= ",Ea.m.REa/Ea.m.DEa))
	print(paste0("Ea-PPE= ",Ea.m.RPPE/Ea.m.DPPE))
print("############################################################")	
	Ps.m.DCtl=mean(Ps[Ps$Concat =="Dcontrol"]$relabund)
	Ps.m.DEa=mean(Ps[Ps$Concat =="DEa"]$relabund)
	Ps.m.DPPE=mean(Ps[Ps$Concat =="DPPE"]$relabund)
	
	Ps.m.RCtl=mean(Ps[Ps$Concat =="Rcontrol"]$relabund)
	Ps.m.REa=mean(Ps[Ps$Concat =="REa"]$relabund)
	Ps.m.RPPE=mean(Ps[Ps$Concat =="RPPE"]$relabund)	
	print(paste0("Ps-Ctl= ",Ps.m.RCtl/Ps.m.DCtl))
	print(paste0("Ps-Ea= ",Ps.m.REa/Ps.m.DEa))
	print(paste0("Ps-PPE= ",Ps.m.RPPE/Ps.m.DPPE))
print("############################################################")		
	Pa.m.DCtl=mean(Pa[Pa$Concat =="Dcontrol"]$relabund)
	Pa.m.DEa=mean(Pa[Pa$Concat =="DEa"]$relabund)
	Pa.m.DPPE=mean(Pa[Pa$Concat =="DPPE"]$relabund)
	
	Pa.m.RCtl=mean(Pa[Pa$Concat =="Rcontrol"]$relabund)
	Pa.m.REa=mean(Pa[Pa$Concat =="REa"]$relabund)
	Pa.m.RPPE=mean(Pa[Pa$Concat =="RPPE"]$relabund)
	print(paste0("Pa-Ctl= ",Pa.m.RCtl/Pa.m.DCtl))
	print(paste0("Pa-Ea= ",Pa.m.REa/Pa.m.DEa))
	print(paste0("Pa-PPE= ",Pa.m.RPPE/Pa.m.DPPE))
print("############################################################")		
	
	
	
	
	Ea.Ctl=data.table(Ea[Ea$Concat == "Rcontrol"]$relabund/Ea.m.DCtl)
	Ea.Ctl$Bac="Ea"
	Ea.Ctl$Treat="Ctl"
	
	Ea.Ea=data.table(Ea[Ea$Concat == "REa"]$relabund/Ea.m.DEa)
	Ea.Ea$Bac="Ea"
	Ea.Ea$Treat="Ea"
	
	Ea.PPE=data.table(Ea[Ea$Concat == "RPPE"]$relabund/Ea.m.DPPE)
	Ea.PPE$Bac="Ea"
	Ea.PPE$Treat="PPE"
	
	Ps.m.DCtl=mean(Ps[Ps$Concat =="Dcontrol"]$relabund)
	Ps.m.DEa=mean(Ps[Ps$Concat =="DEa"]$relabund)
	Ps.m.DPPE=mean(Ps[Ps$Concat =="DPPE"]$relabund)
	
	Ps.Ctl=data.table(Ps[Ps$Concat == "Rcontrol"]$relabund/Ps.m.DCtl)
	Ps.Ctl$Bac="Ps"
	Ps.Ctl$Treat="Ctl"
	
	Ps.Ea=data.table(Ps[Ps$Concat == "REa"]$relabund/Ps.m.DEa)
	Ps.Ea$Bac="Ps"
	Ps.Ea$Treat="Ea"
	
	Ps.PPE=data.table(Ps[Ps$Concat == "RPPE"]$relabund/Ps.m.DPPE)
	Ps.PPE$Bac="Ps"
	Ps.PPE$Treat="PPE"
	
	
	
	Pa.Ctl=data.table(Pa[Pa$Concat == "Rcontrol"]$relabund/Pa.m.DCtl)
	Pa.Ctl$Bac="Pa"
	Pa.Ctl$Treat="Ctl"
	
	Pa.Ea=data.table(Pa[Pa$Concat == "REa"]$relabund/Pa.m.DEa)
	Pa.Ea$Bac="Pa"
	Pa.Ea$Treat="Ea"
	
	Pa.PPE=data.table(Pa[Pa$Concat == "RPPE"]$relabund/Pa.m.DPPE)
	Pa.PPE$Bac="Pa"
	Pa.PPE$Treat="PPE"
	
	DT.m=rbind(Ea.Ea, Pa.Ea, Ps.Ea,Ea.Ctl, Pa.Ctl, Ps.Ctl,Ea.PPE, Pa.PPE, Ps.PPE)
	DT.m$Concat<-paste0(DT.m$Bac,"_",DT.m$Treat)
	
	pp6=ggplot(DT.m, aes(x=Concat, y=V1, color=Bac))
	p6=pp6+ geom_boxplot()+ geom_point(size=2)+
	theme_bw()+ theme_new1 +
	scale_color_manual(values=color.2) +
	ylab("16S rRNA : 16S rDNA ratio")

#	pdf("box_plot_rna_dna.pdf",useDingbats=FALSE)
	print(p6)
#	dev.off()
	
	
	
	
	
	
	
	
	
	
	

