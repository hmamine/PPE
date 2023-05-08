#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "picante",  "PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c("#cc2a36","#166590","#11902d")


theme_new <- theme (
	panel.border = element_rect(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x = element_text(angle=90, vjust=1),
#	axis.text.x=element_blank(),
	axis.ticks.x=element_blank()
)
X=(0)
Y=(1)

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

#physeq.2D=subset_samples( physeq.BT, concat1 != "4D2dpi" )

#trimming OTUs with low occurance
otu.I=otu_table(physeq.BT)
taxa.I=tax_table(physeq.BT)
sd.I=sample_data(physeq.BT)
flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
otu.I.filtered=flt$mat
taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
dummyTaxonomy=c("D__dummy","Ph__","Cl__","Or__","Fa__","Ge__","Sp__","Se__")
TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
NewTaxa=tax_table(TAXA.filt)

physeq.flt=phyloseq(NewOTU, NewTaxa, sd.I)

physeq.flt=tax_glom(physeq.flt, "Family") 
physeq.RA=transform_sample_counts(physeq.flt, function(x) x / sum(x) )

DT_sd=data.table( as( sample_data( physeq.RA ), "data.frame" ), keep.rownames=T )
setnames(DT_sd, "rn", "Samples")
setkey(DT_sd, "Samples")
DT_RA=melt(as(otu_table( physeq.RA ), "matrix" ))
DT_RA=data.table( DT_RA,  keep.rownames=T )
setkey(DT_RA, "Var1", "Var2")

DT_tax=as(tax_table( physeq.RA ), "matrix" )
DT_tax=data.table(DT_tax, keep.rownames=T)
setnames(DT_tax, "rn", "Names")
setkey(DT_tax, "Names")
	
DT_1<-merge(DT_RA, DT_tax, by.x="Var1", by.y="Names")
DT_2<-merge(DT_1, DT_sd, by.x="Var2", by.y="Samples") 
LIST <- c("Microbacteriaceae","Micrococcaceae","Bacillaceae","Beijerinckiaceae"
,"Erwiniaceae","Oxalobacteraceae","Pseudomonadaceae")

DT1=DT_2[DT_2$concat1=="2D2dpi" & DT_2$Family %in% LIST]
	p1=ggplot(DT1, aes(y=value, x=Family, color=Treat))+
	geom_boxplot()+theme_bw()+theme_new+
	scale_color_manual(values=color.1) 

DT2=DT_2[DT_2$concat1=="2D4dpi" & DT_2$Family %in% LIST]
	p2=ggplot(DT2, aes(y=value, x=Family, color=Treat))+
	geom_boxplot()+theme_bw()+theme_new+
	scale_color_manual(values=color.1)

	pdf("relative_abundance.pdf",useDingbats=FALSE)
        gridExtra::grid.arrange(p1,p2, nrow=2, ncol=2)
	dev.off()
#	DT=DT_2[DT_2$Var1=="43fddf1528d4a98928fd8c3a8ac23bfd"]
#	DT=DT_2[DT_2$Var1=="0f18144d308ada95632ab5193d92073f"]
#	DT=DT_2[DT_2$Var1=="2a39f0b75d31d00b3136b2edbe82962e"]
#	DT=DT_2[DT_2$Var1=="43fddf1528d4a98928fd8c3a8ac23bfd"]
#	DT=DT_2[DT_2$Var1=="5648dccee530d68ceb3e4d7d22cf8756"]
#	DT=DT_2[DT_2$Var1=="69a7261dfbf942221779304cb06406fc"]
#	DT=DT_2[DT_2$Var1=="9af3467db68cf6063627304cecd46a65"]
#	DT=DT_2[DT_2$Var1=="b03d2e5ac94c7a0d115e5c542e9315e8"]
#	DT=DT_2[DT_2$Var1=="ca6c2b3b469c08212142e8821c65882a"]
#	DT=DT_2[DT_2$Var1=="cc761daf51f27c423da57f3f1f0ff5cc"]
#	DT=DT_2[DT_2$Var1=="d6efd2da2728fd74ded268122ee05036"]
#	DT=DT_2[DT_2$Var1=="d829bee4984f82ffc2453212157caf96"]
#	DT=DT_2[DT_2$Var1=="dddf3f4996401c070a8b7ad646f2620f"]
#	DT=DT_2[DT_2$Var1=="e053ffbf21bc73ea1948801c1acdd5df"]
#	DT=DT_2[DT_2$Var1=="e6f35e978bd7ee5baf1f0e2b955965f6"]
#	DT=DT_2[DT_2$Var1=="f0517e4077d3468f0ceea9b5f1ab55e7"]
#	p1=ggplot(DT, aes(y=value, x=concat2, color=Treat))+
#	geom_boxplot()+geom_point()+theme_bw()+theme_new
#	print(p1)
	
	
	







