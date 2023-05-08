#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "DESeq2", "venn","readr","magrittr")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
color.list1<-c("#77aaff","#ff77aa","#bcbcbc")
alpha = 0.01	
Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#000000","#2BB065")
LIST<-c("Ea_01648_Transcriptional_regulator_SlyA","Ea_00767_Blue_copper_oxidase_CueO",
"Ea_03025_Regulator_of_ribonuclease_activity_B")
Shape.1=c(21,19)

theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.ticks = element_line(colour="black"),
	panel.background = element_blank(), 
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none",
	)

theme_new2 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	legend.position="none",
	)

theme_new3 <- theme(
	axis.text = element_text(colour="black", size=8),
	text = element_text(size=8),
  	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	panel.background = element_blank(),
	axis.ticks = element_line(colour="black"),
	axis.line = element_line(colour="black"),
	legend.position="none",
	)	


Volcan<-function(mapping){
		p<- ggplot( DT2 , mapping )
		p+geom_jitter( size=1 )+ theme_bw( )+
		geom_hline(yintercept = -log2(0.05), linetype="longdash")+
		geom_vline(xintercept = c(-1.5,0,1.5), linetype="longdash")+
		theme_new3
	}

DT<-vector( "list" )
DT.all<-vector( "list" )
pEa<-vector( "list" )
pPa<-vector( "list" )
pPs<-vector( "list" )
pMA.Ea<-vector( "list" )
pMA.Pa<-vector( "list" )
pMA.Ps<-vector( "list" )
#log2 fold change
	folch=1.5
#Adjusted p-value
	alpha=0.05
	
 KO=read.delim("KO_terms.txt",sep="\t", row.names=1, header=T)
 tmp=read.table("tmp.txt",sep="\t", row.names=1, header=T)
# upload and prepare phyloseq objects***
	tab=read.table("Uni_quant_count.sf", sep="\t", row.names=1, header=T)
	mat=as.matrix(tab)
	OTU=otu_table(mat, taxa_are_rows=T) 
	sd=read.table("Uni_quant_count.sd", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)	
	tx=read.table("Uni_quant.id", sep="\t", row.names=1, header=1)
	tx=as.matrix(tx)
	TX=tax_table(tx)
	physeq=phyloseq(OTU, TX, SD ) 
	
###computing differential expressed genes in Erwinia
	print ("###computing differential expressed genes in Erwinia")
	physeq1=subset_taxa( physeq, Strain == "Ea" )
	LIST1<-list ("PaEa","PsEa","PPE")
	DT<-vector( "list" )
	for( i in LIST1 ) 
	{
		print( i )
		SUBSET=c("Ea", i )
		physeq2=subset_samples( physeq1, Cond %in% SUBSET )
		physeq.i = filter_taxa(physeq2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
		dds.i=phyloseq_to_deseq2(physeq.i, ~Cond)
		dds.i=DESeq(dds.i, fitType='parametric')
		res.i = results(dds.i, cooksCutoff = FALSE, tidy=TRUE)
		rownames(res.i)<-res.i$row
		res.i<-res.i[,-1]
		
		DT_res<-data.table(res.i, keep.rownames=T, key="rn")
		DT.taxa=data.table( as(tax_table(physeq.i),"matrix"), keep.rownames=T, key="rn")
		DT.SD=data.table( as(sample_data(physeq.i),"data.frame"), keep.rownames=T, key="rn")
		DT_i=merge(DT_res, DT.taxa, by.x="rn", by.y="rn")

		DT2=DT_i
		p1<-Volcan( aes( y=-log2(padj), x=log2FoldChange)) +
		scale_y_continuous(limits=c(0,160))+
		scale_x_continuous(limits=c(-6,6), breaks=c(-5.5,-3.5,-1.5,0,1.5,3.5,5.5))
		
#		print(p1)
		pEa[[i]]=p1
	
#		plotCounts(dds.i, gene="Ea_03344_hypothetical_protein", intgroup="Cond")
		DT[[i]]=DT2[DT2$padj < 0.05 & abs(DT2$log2FoldChange) > 1.5,]

		m=as( otu_table( physeq.i ), "matrix" )
		DT.ra=data.table(m, keep.rownames=T)
		DT.ra$mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = colnames(m)]
		DT.ra=DT.ra[,c("rn","mean")]	

		DT1=merge(DT.ra, DT_i, by.x="rn",by.y="rn")
#		DT1$Col<-ifelse(DT1$log2FoldChange < -1.5 & DT1$padj < 0.05, "Dep.", ifelse(DT1$log2FoldChange > 1.5 & DT1$padj < 0.05, "Enr.", "N.S."))

		DT1$Col<-ifelse(DT1$log2FoldChange < -1.5 & DT1$padj < 0.05, "Dep.", ifelse(DT1$log2FoldChange > 1.5 & DT1$padj < 0.05, "Enr.", "N.S."))
		
		DT1$Nam<-ifelse(!DT1$rn %in% LIST, "",DT1$rn)
		
		pp2=ggplot(DT1, aes(x=log2(mean), y=log2FoldChange, color=Col))
		p2=pp2+geom_point(size=2.5)+
		scale_y_continuous(limits=c(-8,4), breaks=c(-7.5,-5.5,-3.5,-1.5,0,1.5,3.5))+
		scale_x_continuous(limits=c(0,25))+
		theme_bw()+theme_new1+geom_text(aes(label=Nam))+
		scale_color_manual(values=color.list1)
		
#		print(p2)
		pMA.Ea[[i]]=p2	
	}

#	pdf("MA_plot_Ea.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(pMA.Ea$PaEa, pMA.Ea$PsEa, pMA.Ea$PPE, nrow=2, ncol=2)
#	dev.off()

	m=as( otu_table( physeq1 ), "matrix" )
	DT_m <- m %>% 
	melt(id.vars = "variable", value.name="relabund")
	DT.m=data.table(DT_m, keep.rownames=F, key="Var2") 
	
	DT.taxa=data.table( as(TX,"matrix"), keep.rownames=T, key="rn")
	DT.SD=data.table( as(SD,"data.frame"), keep.rownames=T, key="rn")
	
	DT1=merge(DT.m, DT.taxa, by.x="Var1",by.y="rn")
	DT1=merge(DT1, DT.SD, by.x="Var2",by.y="rn")

	PaEa_Pos=DT$PaEa[DT$PaEa$log2FoldChange > 0,]$rn
	PaEa_Neg=DT$PaEa[DT$PaEa$log2FoldChange < 0,]$rn
	
	PsEa_Pos=DT$PsEa[DT$PsEa$log2FoldChange > 0,]$rn
	PsEa_Neg=DT$PsEa[DT$PsEa$log2FoldChange < 0,]$rn
	
	PPE_Pos=DT$PPE[DT$PPE$log2FoldChange > 0,]$rn
	PPE_Neg=DT$PPE[DT$PPE$log2FoldChange < 0,]$rn

#	pdf("venn_positive_Ea.pdf",useDingbats=FALSE)
#	V1=venn(list(PaEa_Pos,PsEa_Pos,PPE_Pos), zcolor="style", snames=c("Pa","Ps","Both"))
#	ATTR1<-attr(V1,"intersection")
#	dev.off()

#	pdf("venn_negative_Ea.pdf",useDingbats=FALSE)
#	V2=venn(list(PaEa_Neg,PsEa_Neg,PPE_Neg), zcolor="style", snames=c("Pa","Ps","Both"))
#	ATTR2<-attr(V2,"intersection")
#	dev.off()	


###computing differential expressed genes in Pantoea
	print ("###computing differential expressed genes in Pantoea")
	physeq1=subset_taxa( physeq,  Strain == "Pa")
	LIST2<-list ("PaEa","PP","PPE")
	DT<-vector( "list" )
	for( i in LIST2 ) 
	{
		print( i )
		SUBSET=c("Pa", i )

		physeq2=subset_samples( physeq1, Cond %in% SUBSET )
		physeq.i = filter_taxa(physeq2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
		dds.i=phyloseq_to_deseq2(physeq.i, ~Cond)
		dds.i=DESeq(dds.i, fitType='mean')
		res.i = results(dds.i, cooksCutoff = FALSE, tidy=TRUE)
		rownames(res.i)<-res.i$row
		res.i<-res.i[,-1]

		DT_res<-data.table(res.i, keep.rownames=T, key="rn")
		DT.taxa=data.table( as(tax_table(physeq.i),"matrix"), keep.rownames=T, key="rn")
		DT.SD=data.table( as(sample_data(physeq.i),"data.frame"), keep.rownames=T, key="rn")
		DT_i=merge(DT_res, DT.taxa, by.x="rn", by.y="rn")

		DT2=DT_i
		p1<-Volcan( aes( y=-log2(padj), x=log2FoldChange)) +
		scale_y_continuous(limits=c(0,160))+
		scale_x_continuous(limits=c(-6,6), breaks=c(-5.5,-3.5,-1.5,0,1.5,3.5,5.5))
		
#		print(p1)
		pPa[[i]]=p1
	
#		plotCounts(dds.i, gene="Ea_03344_hypothetical_protein", intgroup="Cond")
		DT[[i]]=DT2[DT2$padj < 0.05 & abs(DT2$log2FoldChange) > 1.5,]

		m=as( otu_table( physeq.i ), "matrix" )
		DT.ra=data.table(m, keep.rownames=T)
		DT.ra$mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = colnames(m)]
		DT.ra=DT.ra[,c("rn","mean")]	

		DT1=merge(DT.ra, DT_i, by.x="rn",by.y="rn")
		DT1$Col<-ifelse(DT1$log2FoldChange < -1.5 & DT1$padj < 0.05, "Dep.", ifelse(DT1$log2FoldChange > 1.5 & DT1$padj < 0.05, "Enr.", "N.S."))
		
		pp2=ggplot(DT1, aes(x=log2(mean), y=log2FoldChange, color=Col))
		p2=pp2+geom_point(size=2.5)+
		theme_bw()+theme_new1+
		scale_y_continuous(limits=c(-6,6), breaks=c(-4.5,-3,-1.5,0,1.5,3,4.5))+
		scale_x_continuous(limits=c(0,25))+
		scale_color_manual(values=color.list1)
		
#		print(p2)
		pMA.Pa[[i]]=p2
	}
	
#	pdf("MA_plots_Pa.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(pMA.Pa$PaEa, pMA.Pa$PP, pMA.Pa$PPE, nrow=2, ncol=2)
#	dev.off()

#	PaEa_Pos=DT$PaEa[DT$PaEa$log2FoldChange > 0,]$rn
#	PaEa_Neg=DT$PaEa[DT$PaEa$log2FoldChange < 0,]$rn
	
#	PP_Pos=DT$PP[DT$PP$log2FoldChange > 0,]$rn
#	PP_Neg=DT$PP[DT$PP$log2FoldChange < 0,]$rn
	
#	PPE_Pos=DT$PPE[DT$PPE$log2FoldChange > 0,]$rn
#	PPE_Neg=DT$PPE[DT$PPE$log2FoldChange < 0,]$rn
	
#	pdf("venn_positive_Pa.pdf",useDingbats=FALSE)
#	venn(list(PaEa_Pos,PP_Pos,PPE_Pos), zcolor="style", snames=c("Ea","Ps","Both"))
#	dev.off()

#	pdf("venn_negative_Pa.pdf",useDingbats=FALSE)
#	venn(list(PaEa_Neg,PP_Neg,PPE_Neg), zcolor="style", snames=c("Ea","Ps","Both"))
#	dev.off()	


###computing differential expressed genes in Pseudomonas
	print ("###computing differential expressed genes in Pseudomonas")
	physeq1=subset_taxa( physeq,  Strain == "Ps")
	LIST3<-list ("PsEa","PP","PPE")
	DT<-vector( "list" )
	for( i in LIST3 ) 
	{
		print( i )
		SUBSET=c("Ps", i )

		physeq2=subset_samples( physeq1, Cond %in% SUBSET )
		physeq.i = filter_taxa(physeq2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
		dds.i=phyloseq_to_deseq2(physeq.i, ~Cond)
		dds.i=DESeq(dds.i, fitType='mean')
		res.i = results(dds.i, cooksCutoff = FALSE, tidy=TRUE)
		rownames(res.i)<-res.i$row
		res.i<-res.i[,-1]
		
		DT_res<-data.table(res.i, keep.rownames=T, key="rn")
		DT.taxa=data.table( as(tax_table(physeq.i),"matrix"), keep.rownames=T, key="rn")
		DT.SD=data.table( as(sample_data(physeq.i),"data.frame"), keep.rownames=T, key="rn")
		DT_i=merge(DT_res, DT.taxa, by.x="rn", by.y="rn")

		DT2=DT_i
		p1<-Volcan( aes( y=-log2(padj), x=log2FoldChange)) +
		scale_y_continuous(limits=c(0,160))+
		scale_x_continuous(limits=c(-6,6), breaks=c(-5.5,-3.5,-1.5,0,1.5,3.5,5.5))
		
#		print(p1)
		pPs[[i]]=p1
	
#		plotCounts(dds.i, gene="Ea_03344_hypothetical_protein", intgroup="Cond")
		DT[[i]]=DT2[DT2$padj < 0.05 & abs(DT2$log2FoldChange) > 1.5,]

		m=as( otu_table( physeq.i ), "matrix" )
		DT.ra=data.table(m, keep.rownames=T)
		DT.ra$mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = colnames(m)]
		DT.ra=DT.ra[,c("rn","mean")]	

		DT1=merge(DT.ra, DT_i, by.x="rn",by.y="rn")
		DT1$Col<-ifelse(DT1$log2FoldChange < -1.5 & DT1$padj < 0.05, "Dep.", ifelse(DT1$log2FoldChange > 1.5 & DT1$padj < 0.05, "Enr.", "N.S."))
		
		pp2=ggplot(DT1, aes(x=log2(mean), y=log2FoldChange, color=Col))
		p2=pp2+geom_point(size=2.5)+
		theme_bw()+theme_new1+
		scale_y_continuous(limits=c(-6,6), breaks=c(-4.5,-3,-1.5,0,1.5,3,4.5))+
		scale_x_continuous(limits=c(0,25))+
		scale_color_manual(values=color.list1)
		
		print(p2)
		pMA.Ps[[i]]=p2
	}
	
#	pdf("MA_plot_Ps.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(pMA.Ps$PsEa, pMA.Ps$PP, pMA.Ps$PPE, nrow=2, ncol=2)
#	dev.off()

#	PsEa_Pos=DT$PsEa[DT$PsEa$log2FoldChange > 0,]$rn
#	PsEa_Neg=DT$PsEa[DT$PsEa$log2FoldChange < 0,]$rn
#	
#	PP_Pos=DT$PP[DT$PP$log2FoldChange > 0,]$rn
#	PP_Neg=DT$PP[DT$PP$log2FoldChange < 0,]$rn
	
#	PPE_Pos=DT$PPE[DT$PPE$log2FoldChange > 0,]$rn
#	PPE_Neg=DT$PPE[DT$PPE$log2FoldChange < 0,]$rn
	
#	pdf("venn_positive_Ps.pdf",useDingbats=FALSE)
#	V1=venn(list(PsEa_Pos,PP_Pos,PPE_Pos), zcolor="style", snames=c("Ea","Pa","Both"))
#	ATTR1<-attr(V1,"intersection")
#	dev.off()

#	pdf("venn_negative_Ps.pdf",useDingbats=FALSE)
#	venn(list(PsEa_Neg,PP_Neg,PPE_Neg), zcolor="style", snames=c("Ea","Pa","Both"))
#	dev.off()	


















