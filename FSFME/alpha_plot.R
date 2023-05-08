#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "picante",  "PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")
color_palette<-c("#cc2a36","#166590","#11902d","#000000")


theme_new <- theme (
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.ticks.x=element_blank(),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
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

physeq= phyloseq(OTU, TAXA, SD,TREE)

physeq.B=subset_samples( physeq, Treat != "Negative" )
physeq.BM=subset_taxa( physeq.B, !Genus %in% c("Mitochondria","Chloroplast") )
physeq.BT=subset_taxa( physeq.BM, Select != "Prune" )


#Rarefication to even depth based on smallest sample size
#	rf=min(sample_sums(physeq.BT))
#	physeq.rrf=rarefy_even_depth(physeq.BT, 300, replace=TRUE, rngseed = 131)

#Ploting species richness	
	p_nat=plot_richness(physeq.BT,"concat2","Treat" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]

	neworder=c("untreat","2D2dpi_H2O","2D2dpi_Ea","2D4dpi_H2O","2D4dpi_Ea","4D2dpi_H2O","4D2dpi_Ea")
	p_nat$data$concat2 <- ordered(p_nat$data$concat2, levels=neworder )
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=concat2, y=value, color=Treat))+
	geom_point(size=4)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)

subp1 =	ggplot(data=p_nat$data[p_nat$data$variable=="Observed",],aes(x=concat2, y=value, color=Treat))+
	geom_boxplot(width=0.5)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) + ylab("Observed ASVs")

subp2 =ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",],aes(x=concat2, y=value, color=Treat))+
	geom_boxplot(width=0.5)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Shannon Index")

otu.mat.rrf=as( otu_table( physeq.BT), "matrix" )
Faith.PD<-data.table(pd(t(otu.mat.rrf), TREE), keep.rownames=T, key="rn")
setnames(Faith.PD, "rn", "samples")						
DT_PD<-merge(subp2$data[,c(-7:-8)], Faith.PD)
DT_PD$variable <- "FaithPD"

subp3 = ggplot(data=DT_PD, aes(x=concat2, y=PD, color=Treat))+
	geom_boxplot(width=0.5)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Phylogenetic diversity")
	
#print(kruskal.test(data=subp3$data[subp3$data$concat1=="2D2dpi",], PD ~ Treat))
#print(kruskal.test(data=subp3$data[subp3$data$concat1=="2D4dpi",], PD ~ Treat))
#print(kruskal.test(data=subp3$data[subp3$data$concat1=="4D2dpi",], PD ~ Treat))

#	pdf("alpha_diversity.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(subp1, subp3, subp3, nrow=2, ncol=2)
#	dev.off()

           
            
          
 









