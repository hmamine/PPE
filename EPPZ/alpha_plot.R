#***require packages for the analysis the analysis
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus","agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")

color_palette<-c("#155832","#ffa500","#00b2b2","#4f7b95")

set.seed(123)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
#	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
#	legend.position="none",
	axis.text.x = element_text(angle=60, vjust=1, hjust=1, size=7, color="black"),
#	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	)
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
}
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
	
#Rarefication to even depth based on smallest sample size
rf=min(sample_sums(physeq.B))
physeq.rrf=rarefy_even_depth(physeq.B, rf, replace=TRUE, rngseed = 131)

#Ploting species richness	
	p_nat=plot_richness(physeq.rrf,"Concat","Treat" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
	
	Ord1<-c("RNA.control","DNA.control","RNA.Strep","DNA.Strep",
	"RNA.Ea","DNA.Ea","RNA.PPE","DNA.PPE")  
	p_nat$data$Concat<-ordered(p_nat$data$Concat, levels=Ord1 )
	 
p1 = p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat, y=value, color=Treat))+
	geom_point(size=.5)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)

subp0 = ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Treat), outlier.shape = 1, outlier.size= 1)+
	scale_fill_manual(values=color_palette) 
	legend <- g_legend(subp0) 
	
subp2 =ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Treat), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette) +
	ylab("Observed ASVs")
#testing significance in observed ASVs
	kruskal.test(data=subp2$data, value ~ Concat)
	kwAllPairsConoverTest(value ~ as.factor(Concat) , data=subp2$data,  p.adjust.method="BH" )


subp1 = ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",], aes(x=Concat, y=value))+
	geom_boxplot(aes(fill=Treat), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_fill_manual(values=color_palette) + 
	ylab("Shannon index")
#testing significance in Shannon index
	kruskal.test(data=subp1$data, value ~ Concat)
	kwAllPairsConoverTest(value ~ as.factor(Concat) , data=subp1$data,  p.adjust.method="BH" )

#	pdf("alpha_diversity.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(subp1, subp2, legend, nrow=2, ncol=2)
#	dev.off()
           

