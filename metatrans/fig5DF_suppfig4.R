#require packages for the analysis the analysis
pkg=c("ggplot2", "data.table", "reshape2","PMCMR","venn","dplyr","packcircles","scales","tidyr","igraph","RColorBrewer")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])


list.color<-c("#b02b76","#ffa500")
list.color1<-c("#ffa500","#6c6b20","#b36200")
list.color2<-c("#0072ff","#398000","#800039")
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
#	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.y=element_blank(),
#	axis.text.y=element_text(size=8, color="black"),
	axis.ticks.y=element_blank(),
	)

theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
#	axis.text.y=element_blank(),
	axis.text.y=element_text(color="black", size=8)
	)


	DEG=read.delim("diff_enri_all.txt", sep="\t", row.names=1, header=1)
	DT.DEG=data.table(DEG, key="Id")

	DT.Ea=DT.DEG[DT.DEG$strain =="Ea"]
	
	Ea.Pos.coPs=DT.Ea[DT.Ea$coCulture=="coPs" & DT.Ea$log2FoldChange>0]$Id
	Ea.Neg.coPs=DT.Ea[DT.Ea$coCulture=="coPs" & DT.Ea$log2FoldChange<0]$Id
	
	Ea.Pos.coPa=DT.Ea[DT.Ea$coCulture=="coPa" & DT.Ea$log2FoldChange>0]$Id
	Ea.Neg.coPa=DT.Ea[DT.Ea$coCulture=="coPa" & DT.Ea$log2FoldChange<0]$Id
	
	Ea.Pos.coPaPs=DT.Ea[DT.Ea$coCulture=="coPaPs" & DT.Ea$log2FoldChange>0]$Id
	Ea.Neg.coPaPs=DT.Ea[DT.Ea$coCulture=="coPaPs" & DT.Ea$log2FoldChange<0]$Id
	
#	V1=venn(list(Ea.Neg.coPa,Ea.Neg.coPs,Ea.Neg.coPaPs), zcolor="style", snames=c("Pa","Ps","PPE"))
#	ATTR1<-attr(V1,"intersection")	
#	V2=venn(list(Ea.Pos.coPa,Ea.Pos.coPs,Ea.Pos.coPaPs), zcolor="style", snames=c("Pa","Ps","PPE"))
#	ATTR2<-attr(V2,"intersection")


	DT.Pa=DT.DEG[DT.DEG$strain =="Pa"]
	
	Pa.Pos.coPs=DT.Pa[DT.Pa$coCulture=="coPs" & DT.Pa$log2FoldChange>0]$Id
	Pa.Neg.coPs=DT.Pa[DT.Pa$coCulture=="coPs" & DT.Pa$log2FoldChange<0]$Id
	
	Pa.Pos.coEa=DT.Pa[DT.Pa$coCulture=="coEa" & DT.Pa$log2FoldChange>0]$Id
	Pa.Neg.coEa=DT.Pa[DT.Pa$coCulture=="coEa" & DT.Pa$log2FoldChange<0]$Id
	
	Pa.Pos.coEaPs=DT.Pa[DT.Pa$coCulture=="coEaPs" & DT.Pa$log2FoldChange>0]$Id
	Pa.Neg.coEaPs=DT.Pa[DT.Pa$coCulture=="coEaPs" & DT.Pa$log2FoldChange<0]$Id
		
#	V3=venn(list(Pa.Neg.coEa,Pa.Neg.coPs,Pa.Neg.coEaPs), zcolor="style", snames=c("Ea","Ps","PPE"))
#	ATTR3<-attr(V3,"intersection")
#	V4=venn(list(Pa.Pos.coEa,Pa.Pos.coPs,Pa.Pos.coEaPs), zcolor="style", snames=c("Ea","Ps","PPE"))
#	ATTR4<-attr(V4,"intersection")
	

	DT.Ps=DT.DEG[DT.DEG$strain =="Ps"]
	
	Ps.Pos.coPa=DT.Ps[DT.Ps$coCulture=="coPa" & DT.Ps$log2FoldChange>0]$Id
	Ps.Neg.coPa=DT.Ps[DT.Ps$coCulture=="coPa" & DT.Ps$log2FoldChange<0]$Id
	
	Ps.Pos.coEa=DT.Ps[DT.Ps$coCulture=="coEa" & DT.Ps$log2FoldChange>0]$Id
	Ps.Neg.coEa=DT.Ps[DT.Ps$coCulture=="coEa" & DT.Ps$log2FoldChange<0]$Id
	
	Ps.Pos.coEaPa=DT.Ps[DT.Ps$coCulture=="coEaPa" & DT.Ps$log2FoldChange>0]$Id
	Ps.Neg.coEaPa=DT.Ps[DT.Ps$coCulture=="coEaPa" & DT.Ps$log2FoldChange<0]$Id
	
#	V5=venn(list(Ps.Neg.coEa,Ps.Neg.coPa,Ps.Neg.coEaPa), zcolor="style", snames=c("Ea","Pa","PPE"))
#	ATTR5<-attr(V5,"intersection")
#	V6=venn(list(Ps.Pos.coEa,Ps.Pos.coPa,Ps.Pos.coEaPa), zcolor="style", snames=c("Ea","Pa","PPE"))
#	ATTR6<-attr(V6,"intersection")	
	
#######################################################################################

	DEa=read.table("DE_Ea_kodown.txt", sep="\t", row.names=1, header=1)
	DEa1=DEa[c(2,4)]
	
	DPa=read.table("DE_Pa_kodown.txt", sep="\t", row.names=1, header=1)
	DPa1=DPa[c(2,4)]
	
	DPs=read.table("DE_Ps_kodown.txt", sep="\t", row.names=1, header=1)
	DPs1=DPs[c(2,4)]
	
	DPaPs=merge(DPa1,DPs1, by="ko_term", all=TRUE)
	DPPE=merge(DPaPs, DEa1, by="ko_term", all=TRUE)
	DPPE[is.na(DPPE)]<-"0"
	
	name_Ea=data.table(DEa[(1:3)],keep.rownames=T, key="rn")
	setkey(name_Ea, "ko_term")
	name_Pa=data.table(DPa[(1:3)],keep.rownames=T, key="rn")
	setkey(name_Pa, "ko_term")
	name_Ps=data.table(DPs[(1:3)],keep.rownames=T, key="rn")
	setkey(name_Ps, "ko_term")
	
	names=rbind(name_Ea, name_Pa, name_Ps)
	DT.names=names %>% distinct(ko_term, .keep_all = TRUE)

	name_PaPs=merge(name_Pa,name_Ps, by.x="ko_term", by.y="ko_term")
	
	UEa=read.table("DE_Ea_koup.txt", sep="\t", row.names=1, header=1)
	UEa1=UEa[c(2,4)]
	
	UPa=read.table("DE_Pa_koup.txt", sep="\t", row.names=1, header=1)
	UPa1=UPa[c(2,4)]
	
	UPs=read.table("DE_Ps_koup.txt", sep="\t", row.names=1, header=1)
	UPs1=UPs[c(2,4)]
	
	UPaPs=merge(UPa1,UPs1, by="ko_term", all=TRUE)
	UPPE=merge(UPaPs, UEa1, by="ko_term", all=TRUE)
	UPPE[is.na(UPPE)]<-"0"
	
	UPPE=merge(UPaPs, UEa1, by="ko_term", all=TRUE)
	UPPE[is.na(UPPE)]<-"0"
	
	PPE=merge(UPPE, DPPE, by="ko_term", all=TRUE)
	PPE[is.na(PPE)]<-"0"
	rownames(PPE)=PPE$ko_term

	PPE$Down.Ea=as.numeric(PPE$Down.Ea)*-1
 	PPE$Down.Pa=as.numeric(PPE$Down.Pa)*-1
 	PPE$Down.Ps=as.numeric(PPE$Down.Ps)*-1
 	
 	PPE$Up.Ea=as.numeric(PPE$Up.Ea)
 	PPE$Up.Pa=as.numeric(PPE$Up.Pa)
 	PPE$Up.Ps=as.numeric(PPE$Up.Ps)
 	
 	PPE$UpSum=rowSums(PPE[2:4])
 	PPE$DownSum=rowSums(PPE[5:7])
 	
 	m.PPE=melt(PPE[1:7])
	DT.PPE=data.table(m.PPE, keep.rownames=T, key="rn")
	DT.PPE$Iso=ifelse(DT.PPE$variable %in% c("Up.Ea","Down.Ea"),"Ea", ifelse(DT.PPE$variable %in% c("Up.Pa","Down.Pa"), "Pa","Ps"))
	DT.PPE=DT.PPE[DT.PPE$value != "0"]
	DT.PPE$FC<-ifelse(DT.PPE$value > 0, "Enr","Dep")
	DT.merge=merge(DT.PPE, DT.names, by.x="ko_term",by.y="ko_term")
	
	DT.Ea.Enr=DT.merge[DT.merge$FC == "Enr" & DT.merge$ko_term !="ND" & DT.merge$Iso == "Ea"]
	DT.Ea.Enr=DT.Ea.Enr[order(DT.Ea.Enr$value),]
	pp2=ggplot(DT.Ea.Enr,aes(y=Descript, x=value, fill=prot_fam))
	pp2$data$Descript <- ordered(pp2$data$Descript, levels=rev(DT.Ea.Enr$Descript) )
	p2=pp2+geom_col() + theme_bw()+ theme_new + scale_fill_manual(values=list.color2)
	
	DT.Ea.Dep=DT.merge[DT.merge$FC == "Dep" & DT.merge$ko_term !="ND" & DT.merge$Iso == "Ea"]
	DT.Ea.Dep=DT.Ea.Dep[order(DT.Ea.Dep$value),]
	pp3=ggplot(DT.Ea.Dep,aes(y=Descript, x=abs(value), fill=prot_fam))
	pp3$data$Descript <- ordered(pp3$data$Descript, levels=DT.Ea.Dep$Descript )
	p3=pp3+geom_col() + theme_bw()+ theme_new + scale_fill_manual(values=list.color2)
	
#	gridExtra::grid.arrange(p2,p3,nrow=2, ncol=2)
	
	DT.Pa.Enr=DT.merge[DT.merge$FC == "Enr" & DT.merge$ko_term !="ND" & DT.merge$Iso == "Pa"]
	DT.Pa.Enr=DT.Pa.Enr[order(DT.Pa.Enr$value),]
	pp4=ggplot(DT.Pa.Enr,aes(y=Descript, x=value, fill=prot_fam))
	pp4$data$Descript <- ordered(pp4$data$Descript, levels=rev(DT.Pa.Enr$Descript) )
	p4=pp4+geom_col() + theme_bw()+ theme_new + scale_fill_manual(values=list.color2)
	
	DT.Pa.Dep=DT.merge[DT.merge$FC == "Dep" & DT.merge$ko_term !="ND" & DT.merge$Iso == "Pa"]
	DT.Pa.Dep=DT.Pa.Dep[order(DT.Pa.Dep$value),]
	pp5=ggplot(DT.Pa.Dep,aes(y=Descript, x=abs(value), fill=prot_fam))
	pp5$data$Descript <- ordered(pp5$data$Descript, levels=DT.Pa.Dep$Descript )
	p5=pp5+geom_col() + theme_bw()+ theme_new + scale_fill_manual(values=list.color2)
	
#	gridExtra::grid.arrange(p4,p5,nrow=2, ncol=2)
		
	PPE=PPE[order(PPE$UpSum),]
	DT.Enr=DT.merge[DT.merge$FC == "Enr" & DT.merge$ko_term !="ND"]

	ppE=ggplot(DT.Enr,aes(y=ko_term, x=value, fill=Iso))
	ppE$data$ko_term <- ordered(ppE$data$ko_term, levels=rev(PPE$ko_term) )
	pE=ppE+geom_col() +
	theme_bw()+ theme_new + 
	scale_fill_manual(values=list.color1)+
	geom_text(aes(label=Descript),color="black", size=1)

	PPE=PPE[order(PPE$DownSum),]
	DT.Dep=DT.merge[DT.merge$FC == "Dep" & DT.merge$ko_term !="ND"]

	ppD=ggplot(DT.Dep,aes(y=ko_term, x=abs(value), fill=Iso))
	ppD$data$ko_term <- ordered(ppD$data$ko_term, levels=PPE$ko_term )
	pD=ppD+geom_col() +
	scale_fill_manual(values=list.color1)+
	geom_text(aes(label=Descript), color="black",size=1)+
	theme_bw()+ theme_new  
	
	tab1=read.table("tab.txt", sep="\t", row.names=1, header=1)
	DT.tab1=data.table(tab1, keep.rownames=T, key="rn")

	DT.Ea=DT.tab1[DT.tab1$coCulture == "coPaPs" & DT.tab1$KO_term != "ND"]
	PEa = DT.Ea[DT.Ea$log2FoldChange > 0]$KO_term
	NEa = DT.Ea[DT.Ea$log2FoldChange < 0]$KO_term
		
	DT.Pa=DT.tab1[DT.tab1$coCulture == "coEaPs" & DT.tab1$KO_term != "ND"]
	PPa = DT.Pa[DT.Pa$log2FoldChange > 0]$KO_term
	NPa = DT.Pa[DT.Pa$log2FoldChange < 0]$KO_term
	
	DT.Ps=DT.tab1[DT.tab1$coCulture == "coEaPa" & DT.tab1$KO_term != "ND"]
	PPs = DT.Ps[DT.Ps$log2FoldChange > 0]$KO_term
	NPs = DT.Ps[DT.Ps$log2FoldChange < 0]$KO_term
	
	V=venn(list(PEa,PPa,PPs), zcolor="style", snames=c("Ea","Pa","Ps"))
	AT<-attr(V,"intersection")

	DT.PPE=DT.tab1[DT.tab1$coCulture %in% c("coPaPs","coEaPa","coEaPs") & DT.tab1$KO_term != "ND"]
	
	ABC=c("K23508","K02038", "K02046", "K15554", "K02036", "K23163", "K02047","K10003", 
	"K02045", "K02048", "K02020", "K10111", "K07122", "K02067", "K07323", "K09817", 
	"K02042","K15580", "K10537", "K10539", "K10538", "K23058", "K12368", "K23055",
	"K06160")

	DT.ABC=DT.PPE[DT.PPE$KO_term %in% ABC]
	
	pp.ABC=ggplot(DT.ABC, aes(y=KO_term, x=coCulture))
	
	pp.ABC$data$KO_term <- ordered(pp.ABC$data$KO_term, levels=ABC )
	
	p.ABC=pp.ABC+geom_tile(aes(fill=log2FoldChange, color=log2FoldChange))+
	scale_fill_gradient2(low = "darkblue", high = "darkred", space = "Lab")+
	scale_color_gradient2(low = "darkblue", high = "darkred",  space = "Lab")+
	theme_minimal() + theme_new1
	
#	pdf("plot_ABC_transportors.pdf",useDingbats=FALSE)
#	print(p.ABC)
#	dev.off()	
	
	PDS=c("K18921","K18843","K07319","K07154","K19048","K06218","K07334", 
	"K18923","K18918","K21498","K07339","K19166","K04487","K07062","K03531")
	
	DT.PDS=DT.PPE[DT.PPE$KO_term %in% PDS]
	
	pp.PDS=ggplot(DT.PDS, aes(y=KO_term, x=coCulture))
	
	pp.PDS$data$KO_term <- ordered(pp.PDS$data$KO_term, levels=PDS )
	
	p.PDS=pp.PDS+geom_tile(aes(fill=log2FoldChange, color=log2FoldChange))+
	scale_fill_gradient2(low = "darkblue", high = "darkred", space = "Lab")+
	scale_color_gradient2(low = "darkblue", high = "darkred",  space = "Lab")+
	theme_minimal() + theme_new1

#	pdf("plot_Prot_def_sys.pdf",useDingbats=FALSE)
#	print(p.PDS)
#	dev.off()
	
	SeSy=c("K03210","K03072","K02412","K02409","K02411","K03075","K12072","K02680",
	"K12064","K12063","K12056","K12065","K12067","K12059","K12066","K12068",
	"K02679","K07345","K02282","K12510","K02283","K11890","K12511","K11900","K11895",
	"K11902","K03071","K11901","K11897")
	
	DT.SeSy=DT.PPE[DT.PPE$KO_term %in% SeSy]
	
	pp.SeSy=ggplot(DT.SeSy, aes(y=KO_term, x=coCulture))
	
	pp.SeSy$data$KO_term <- ordered(pp.SeSy$data$KO_term, levels=SeSy )
	
	p.SeSy=pp.SeSy+geom_tile(aes(fill=log2FoldChange, color=log2FoldChange))+
	scale_fill_gradient2(low = "darkblue", high = "darkred", space = "Lab")+
	scale_color_gradient2(low = "darkblue", high = "darkred",  space = "Lab")+
	theme_minimal() + theme_new1

#	pdf("plot_sec_sys.pdf",useDingbats=FALSE)
#	print(p.SeSy)
#	dev.off()

	gridExtra::grid.arrange(pD,pE,nrow=2, ncol=2)
	dev.new()
	gridExtra::grid.arrange(p.ABC,p.PDS,p.SeSy,nrow=2, ncol=2)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
