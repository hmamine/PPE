#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools", "metagenomeSeq","PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
#color.list1<-c("#b02b76","#ffa500","#00b2b2","#155832")
color.list1<-c("#e6194B","#f58231","#4363d8","#2a7d34","#4f7b95","#b36200","#911eb4")
color.list2<-c("#000000","#e6194B","#f58231","#b36200")

theme_new <- theme (
	panel.border = element_rect(),
#	legend.position="none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	)
my.theme <- theme(axis.text = element_text(colour="black", size=10),
	text = element_text(size=10),
        title = element_text(size=10, face="bold", vjust=2),
	panel.background = element_rect(fill = 'gray99',colour = "black", size=0.5),
	axis.title.x=  element_text(vjust=-0.45),
	axis.title.y = element_text(vjust=1.2),
	axis.ticks = element_line(colour="black"),
	axis.line = element_line(),
	panel.grid.major = element_line(colour = "gray40", linetype="dotted"),
	panel.grid.minor = element_line(colour = "gray40", linetype="dashed"),
	legend.justification=c(1,1),
#	legend.position="none",
	legend.title=element_blank(),
	legend.text = element_text(size = 10))	
	
meth2=("NMDS")
meth1=("PCoA")
meth3=("CAP")
dist1=("bray")

# upload and prepare phyloseq objects***
	tab=read.table( "Uni_quant.sf", sep="\t", row.names=1, header=T)
	mat=as.matrix(tab)
	OTU=otu_table(mat, taxa_are_rows=T) 
	sd=read.table("Uni_quant.sd", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)	
	tx=read.table("Uni_quant.id", sep="\t", row.names=1, header=1)
	tx=as.matrix(tx)
	TX=tax_table(tx)
	physeq=phyloseq(OTU, TX, SD ) 
#computing Bray-Curtis distances
	dist_BC=distance( physeq, dist1 )
#computing unconstrained PCoA
	pcoa_BC=ordinate( physeq, meth1 , dist_BC )
	pBC<-plot_ordination( physeq, pcoa_BC, color="Cond")
	pBC$layers<-pBC$layers[-1]
	p1=pBC+geom_point(size=4)+theme_bw()+theme_new+ 
	scale_color_manual(values=color.list1) 

	res=read.table( "Uni_quant.sf", sep="\t", row.names=1, header=T)
	tab=read.table( "Uni_quant.sd", sep="\t", row.names=1, header=T)
	DT_sd<-data.table(tab, keep.rownames=T, key="rn")
	pca_Ea <- prcomp(t(res))
	res_Ea <- summary(pca_Ea)
	scree.Ea <- as.data.frame(res_Ea$importance)
	score.Ea <- as.data.frame(res_Ea$x)
	loadings.Ea <- as.data.frame(res_Ea$rotation)
	data_Ea<-score.Ea[, c(1:4)]
	var1_Ea <-round(scree.Ea[2,1:4] * 100, 1)
	DT_Ea<-data.table(data_Ea, keep.rownames=T, key="rn")
	DT_Ea=DT_Ea[DT_sd]

	pp2<-ggplot(DT_Ea, aes(PC1, PC2))
	p2 <- pp2 + geom_point(aes(colour=Cond), size=4) +
	scale_color_manual(values=color.list1)+
	xlab(paste0("PC1 [",var1_Ea$PC1,"%]")) + 
	ylab(paste0("PC2 [",var1_Ea$PC2,"%]")) + 
	theme_bw() + theme_new

#	nmds_BC=ordinate( physeq, meth2 , dist_BC )
#	plot_ordination( physeq, nmds_BC, color="Cond")

	p3<-plot_ordination(physeq, pcoa_BC, type="biplot",  color="Strain") +
	geom_point(size=1, alpha=.5)+theme_bw()+theme_new +
	scale_color_manual(values=color.list2) 

	p4=ggplot(data=p3$data[p3$data$id.type=="Taxa",], aes(x=Axis.1, y=Axis.2, color="Strain")) + 
	geom_point(size=1, alpha=.7)+
	theme_bw()+theme_new
		
#	pdf("PC_PCO_plots_unidex.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
#	dev.off()

	
#	pdf("Biplot_unidex.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1,nrow=2, ncol=2)
#	dev.off()














