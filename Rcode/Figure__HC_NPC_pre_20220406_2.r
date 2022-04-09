###
  cd /public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2
  
  conda activate r_env
  R
###Library
 ###library
 library(dplyr)
 library(ggplot2)
 require(cowplot)
 library(ggsignif)
 library(ade4)
 library(RColorBrewer)
 library(vegan)
 library(psych)
 library(reshape2)
 library(tidyverse)
 library(VennDiagram)
 library(venn)
 library(UpSetR)
 library(plyr)
 library(magrittr)
 library(randomForest) 
 library(pROC)
 library(readxl)
 library(vegan)
 library(corrplot)
 library(stringr)
 library(ROCR)
 library(data.table)
 library(caret)
 library(lars) ###lasso回归
 library(glmnet) ###lasso回归
 library(circlize)
 library(pheatmap)
 library(ggnewscale)
 library(ggstar)
 library(tidytree)
 library(doBy) 
 library(igraph)#用于共现性网络的绘制
 library(Hmisc)#进行共现性网络构建前的OTU两两相关性计算
 library(GOplot)
 library(ggrepel)
 library(openxlsx)
 library(ggpubr)
 library(MASS) # to access Animals data sets
 library(scales)
 library(ggExtra)
 library(factoextra)
 library(FactoMineR)
 library(survival)
 library(survminer)
 library(indicspecies)#用于基于indicator species的组间差异OTU识别
 library(sciplot)#用于基于模块丰度的绘图
 library(ggpmisc)
 library(OptimalCutpoints)
 library(tableone)
 library(picante)
 library(ape)
 library(pROC)
 library(randomForest)
 library(ggsci)
 library(scatterplot3d)
 library(lavaan)
 library(haven)
 library(partykit)
 library(Boruta)
 library(rpart)
 library(Gmisc) ###显示cluster的转化
 library(phyloseq)
 library(microbiome)
 library(ggtree)
 library(ggtreeExtra)
 library(semPlot)
 library("survival")
 library("survminer")
 library(ggalluvial)
 library(ggtreeExtra)
 library(treeio)
 library(tidytree)
 library(ggridges)
 library(bookdown)
 library(ggstar)
 library(Cairo)
 library(aplot)
 library(patchwork)
 library(ggnewscale)
 library(knitr)
 library(ggpp) 
 library(tibble)
 library(tidyr)
 library(dplyr)
 library(ggtree)
 library(ggimage) ###未成功
 library(kableExtra) ###未成功 依赖stringi
 library(MicrobiotaProcess) ###未成功依赖stringi conda install r-stringi 进行安装
 library(ggpattern)  ###未成功
 library(ape) ###将tree转换成data frame
 library(DOSE)
 library(org.Hs.eg.db)
 library(topGO)
 library(clusterProfiler)
 library(pathview)
 library(MicrobiomeProfiler)
 library(ggExtra)
 
 library(WGCNA)
 library(MLeval) ###评价预测模型 --20211003
 library(caret)
 library(xgboost) 
 library(DALEX) ###模型解释包，搭配caret使用
 library(ggrepel)
 Tax4Fun
 library(partykit)
 library(party)
##inputdata
 ###Input
  #读入数据，整理数据至可用的格式
  NPC311otu<-read.table("NPC311otu-table.txt",header=T,check.names=FALSE)
  NPC311tax<-read.table("NPC311taxonomy.txt",header=T,fill=T,check.names=FALSE)
  NPC311meta <-read.table("NPC311meta.txt",header=T,fill=T,check.names=FALSE)
  NPC311otu.tax <-merge(NPC311tax,NPC311otu)
  head(NPC311otu.tax)
  dim(NPC311otu.tax)
  ###去掉未注释到属的序列
  NPC311otu.tax1<-filter(NPC311otu.tax, Genus != "")
  NPC311otu.tax2<-filter(NPC311otu.tax1, Genus != "g__")
  FeatureID<-NPC311otu.tax2$FeatureID
  FeatureID1<-as.data.frame(FeatureID)
  write.table(FeatureID1,"FeatureID1.txt",sep="\t",quote=F)

  NPC311otu.tax3<-NPC311otu.tax2
  NPC311otu.tax3_Genus<-subset(NPC311otu.tax3,select=-c(Kingdom,Phylum,Class,Order,Family,Species,ControlM1,ControlN1,ControlP1,NPC002CN)) #保留genus水平
  dim(NPC311otu.tax3_Genus)
  NPC311otu.tax3_Genus$sum<-rowSums(NPC311otu.tax3_Genus[,3:ncol(NPC311otu.tax3_Genus)])
  NPC311otu.tax3_Genus1<-NPC311otu.tax3_Genus[!duplicated(NPC311otu.tax3_Genus$Genus), ]
  NPC311otu.tax3_Genus2 <- data.frame(NPC311otu.tax3_Genus1[-1,],check.names=FALSE) #去掉空白行
  NPC311otu.tax3_Genus2 <- data.frame(NPC311otu.tax3_Genus2[-1,],check.names=FALSE) #去掉g__行
  dim(NPC311otu.tax3_Genus2) #661个属
  rownames(NPC311otu.tax3_Genus2)<-NPC311otu.tax3_Genus2$Genus
  NPC311otu.tax3_Genus2t<-as.data.frame(t(NPC311otu.tax3_Genus2)) 
  ####Fecal中的属水平筛选
  HCFandFBTgenus.sum.res1<-as.data.frame(HCFandFBTgenus.sum.res)
  HCFandFBTgenus.sum.res1$mean<-rowSums(HCFandFBTgenus.sum.res1)/ncol(HCFandFBTgenus.sum.res1)
  HCFandFBTgenus.sum.res1_filter0.00005<-filter(HCFandFBTgenus.sum.res1,mean>0.00005)
  dim(HCFandFBTgenus.sum.res1_filter0.00005)
  write.table(HCFandFBTgenus.sum.res1_filter0.00005,"HCFandFBTgenus.sum.res1_filter0.00005.txt",sep="\t",quote=F)

 ####Saliva中的属水平筛选
  HCSandSBTgenus.sum.res1<-as.data.frame(HCSandSBTgenus.sum.res)
  HCSandSBTgenus.sum.res1$mean<-rowSums(HCSandSBTgenus.sum.res1)/ncol(HCSandSBTgenus.sum.res1)
  HCSandSBTgenus.sum.res1_filter0.00005<-filter(HCSandSBTgenus.sum.res1,mean>0.00005)
  dim(HCSandSBTgenus.sum.res1_filter0.00005)
  write.table(HCSandSBTgenus.sum.res1_filter0.00005,"HCSandSBTgenus.sum.res1_filter0.00005.txt",sep="\t",quote=F)
  
 ####NPS中的属水平筛选
  HCNandNBTgenus.sum.res1<-as.data.frame(HCNandNBTgenus.sum.res)
  HCNandNBTgenus.sum.res1$mean<-rowSums(HCNandNBTgenus.sum.res1)/ncol(HCNandNBTgenus.sum.res1)
  HCNandNBTgenus.sum.res1_filter0.00005<-filter(HCNandNBTgenus.sum.res1,mean>0.00005)
  dim(HCNandNBTgenus.sum.res1_filter0.00005)
  write.table(HCNandNBTgenus.sum.res1_filter0.00005,"HCNandNBTgenus.sum.res1_filter0.00005.txt",sep="\t",quote=F)
 ##Feces中的筛选
  HCFandFBTgenus.sum.res1_filter0.001<-filter(HCFandFBTgenus.sum.res1,mean>0.001)
  dim(HCFandFBTgenus.sum.res1_filter0.001) ### 46个属
  write.table(HCFandFBTgenus.sum.res1_filter0.001,"HCFandFBTgenus.sum.res1_filter0.001.txt",sep="\t",quote=F)
 ##Saliva中的筛选
  HCSandSBTgenus.sum.res1_filter0.001<-filter(HCSandSBTgenus.sum.res1,mean>0.001)
  dim(HCSandSBTgenus.sum.res1_filter0.001) ###43个属
  write.table(HCSandSBTgenus.sum.res1_filter0.001,"HCSandSBTgenus.sum.res1_filter0.001.txt",sep="\t",quote=F)
 ###Saliva中的筛选
  HCNandNBTgenus.sum.res1_filter0.001<-filter(HCNandNBTgenus.sum.res1,mean>0.001)
  dim(HCNandNBTgenus.sum.res1_filter0.001) ###69个属
  write.table(HCNandNBTgenus.sum.res1_filter0.001,"HCNandNBTgenus.sum.res1_filter0.001.txt",sep="\t",quote=F)

  ###Fecal_Saliva_NPS中属水平的并集
  Fecal_Saliva_Genus0.001<-union(rownames(HCFandFBTgenus.sum.res1_filter0.001),rownames(HCSandSBTgenus.sum.res1_filter0.001))
  Fecal_Saliva_NPS_Genus0.001<-union(Fecal_Saliva_Genus0.001,rownames(HCNandNBTgenus.sum.res1_filter0.001))
  Fecal_Saliva_NPS_Genus0.001<-gsub("-","_",Fecal_Saliva_NPS_Genus0.001)
  Fecal_Saliva_NPS_Genus0.001_1<-as.data.frame(Fecal_Saliva_NPS_Genus0.001)
  Fecal_Saliva_NPS_Genus0.001_1$Genus<-Fecal_Saliva_NPS_Genus0.001_1$Fecal_Saliva_NPS_Genus0.001
  NPC311otu.tax3_Genus1$Genus<-gsub("-","_",NPC311otu.tax3_Genus1$Genus)
  Fecal_Saliva_NPS_Genus_FeatureID0.001<-merge(Fecal_Saliva_NPS_Genus0.001_1,NPC311otu.tax3_Genus1,by="Genus")
  Fecal_Saliva_NPS_Genus_FeatureID0.001_1<-subset(Fecal_Saliva_NPS_Genus_FeatureID0.001,select=c(FeatureID,Genus))
  write.table(Fecal_Saliva_NPS_Genus_FeatureID0.001_1,"Fecal_Saliva_NPS_Genus_FeatureID0.001_1.txt",sep="\t",quote=F)
  NPC311otu.tax3_Genus_Phylum<-subset(NPC311otu.tax3,select=c(Genus,Phylum))
  NPC311otu.tax3_Genus_Phylum1<-NPC311otu.tax3_Genus_Phylum[!duplicated(NPC311otu.tax3_Genus_Phylum$Genus), ]
  NPC311otu.tax3_Genus_Phylum2 <- data.frame(NPC311otu.tax3_Genus_Phylum1[-1,],check.names=FALSE) #去掉空白行
  NPC311otu.tax3_Genus_Phylum2 <- data.frame(NPC311otu.tax3_Genus_Phylum2[-1,],check.names=FALSE) #去掉g__行
  write.table(NPC311otu.tax3_Genus_Phylum2,"NPC311otu.tax3_Genus_Phylum2.txt",sep="\t",quote=F)
  

  HCFandFBTgenus.sum.res1_filter0.001_2<-HCFandFBTgenus.sum.res1_filter0.001
  FecesandSaliva_genus<-union(rownames(HCFandFBTgenus.sum.res1_filter0.001_2),rownames(HCSandSBTgenus.sum.res1_filter0.001_2))
  NPS_Feces_Saliva_genus<-union(rownames(HCFandFBTgenus.sum.res1_filter0.001_2),FecesandSaliva_genus)
  FecesandSaliva_sampleID<-union(colnames(HCFandFBTgenus.sum.res1_filter0.001_2),colnames(HCSandSBTgenus.sum.res1_filter0.001_2))
  NPS_Feces_Saliva_sampleID<-union(colnames(HCNandNBTgenus.sum.res1_filter0.001_2),FecesandSaliva_sampleID)
  
  NPC311genus2_rownames<-NPC311genus2
  rownames(NPC311genus2_rownames)<-NPC311genus2_rownames$Genus
  NPC311genus2_rownames1<-subset(NPC311genus2_rownames,select=-c(Genus))
  NPC311genus2_rownames1_sampleID<-subset(NPC311genus2_rownames1,select=c(as.character(NPS_Feces_Saliva_sampleID)))
  NPC311genus2_rownames1t<-as.data.frame(t(NPC311genus2_rownames1_sampleID))
  NPC311genus2_rownames1t_filter0.001<-subset(NPC311genus2_rownames1t,select=c(as.character(NPS_Feces_Saliva_genus)))
  NPS_Feces_Saliva_filter0.001<-NPC311genus2_rownames1t_filter0.001

   NPS_Feces_Saliva_filter0.001.dist <- vegdist(NPS_Feces_Saliva_filter0.001,method="bray") #求50Feces个样本的dist，method=“bray”
   #做PcoA图
   NPS_Feces_Saliva_filter0.001.res=pcoa(NPS_Feces_Saliva_filter0.001.dist,correction="none")
   head(NPS_Feces_Saliva_filter0.001.res$value)
   head(NPS_Feces_Saliva_filter0.001.res$vectors)
   NPS_Feces_Saliva_meta<-rbind(HCNandNBTmeta,HCFandFBTmeta,HCSandSBTmeta)

   NPS_Feces_Saliva_meta$sampleType<-c(rep("NPS",79),rep("Feces",50),rep("Saliva",66))
   NPS_Feces_Saliva_meta$Group<-NPS_Feces_Saliva_meta$Treatmenttype
   NPS_Feces_Saliva_meta$Group<-gsub("HCN","HC",NPS_Feces_Saliva_meta$Group)
   NPS_Feces_Saliva_meta$Group<-gsub("HCS","HC",NPS_Feces_Saliva_meta$Group)
   NPS_Feces_Saliva_meta$Group<-gsub("HCF","HC",NPS_Feces_Saliva_meta$Group)
   NPS_Feces_Saliva_meta$Group<-gsub("NBT","NPC",NPS_Feces_Saliva_meta$Group)
   NPS_Feces_Saliva_meta$Group<-gsub("FBT","NPC",NPS_Feces_Saliva_meta$Group)
   NPS_Feces_Saliva_meta$Group<-gsub("SBT","NPC",NPS_Feces_Saliva_meta$Group)

   NPS_Feces_Saliva_filter0.001.ano1 <- adonis(formula = NPS_Feces_Saliva_filter0.001.dist ~ group, data = NPS_Feces_Saliva_meta, permutations = 999)
   
   NPS_Feces_Saliva_filter0.001.P1=as.data.frame(NPS_Feces_Saliva_filter0.001.ano1$aov.tab)["group", "Pr(>F)"]
   NPS_Feces_Saliva_filter0.001.R1= as.data.frame(NPS_Feces_Saliva_filter0.001.ano1$aov.tab)["group", "R2"] ###HCFandFBTR2=0.03524
   NPS_Feces_Saliva_filter0.001.data= NPS_Feces_Saliva_filter0.001.res$vectors[,1:2]
   NPS_Feces_Saliva_filter0.001.data=data.frame(NPS_Feces_Saliva_filter0.001.data)
   NPS_Feces_Saliva_filter0.001.data1<-NPS_Feces_Saliva_filter0.001.data[,-c(3)]

   NPS_Feces_Saliva_filter0.001.data1<-tibble::rownames_to_column(NPS_Feces_Saliva_filter0.001.data1,"sampleID")
  
   colnames(NPS_Feces_Saliva_filter0.001.data)=c('x','y')
   NPS_Feces_Saliva_filter0.001.eig=as.numeric(NPS_Feces_Saliva_filter0.001.res$value[,1])
   NPS_Feces_Saliva_filter0.001.type= NPS_Feces_Saliva_meta$group
   NPS_Feces_Saliva_filter0.001.data2=merge(NPS_Feces_Saliva_filter0.001.data1, NPS_Feces_Saliva_meta,by="sampleID")

   NPS_Feces_Saliva_filter0.001.data2$sampleType<-factor(NPS_Feces_Saliva_filter0.001.data2$sampleType,levels = c("NPS","Saliva","Feces"))
###Figure1
 cd /public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2
  conda activate r_env
  R


  ##Figure1A
  pdf("Figure1A_HCFandFBT_Genus0.001.PCoA_20211229.pdf",width=4,height=3)
   ggplot(NPS_Feces_Saliva_filter0.001.data2,aes(x=x,y=y,color=sampleType))+scale_color_manual(values=c('Feces'="#6495ED","Saliva"="#FFA500","NPS"="#228B22"))+geom_point(alpha=.7,size=2,aes(shape = Group))+scale_shape_manual(values=c(0,16))+labs(x=paste("PCoA1(",format(100*NPS_Feces_Saliva_filter0.001.eig[1]/sum(NPS_Feces_Saliva_filter0.001.eig),digits=4),"%)",sep=""),y=paste("PCoA2(",format(100*NPS_Feces_Saliva_filter0.001.eig[2]/sum(NPS_Feces_Saliva_filter0.001.eig),digits=4),"%)",sep=""),title="")+theme_classic()+stat_ellipse(type="t",linetype=2)+annotate("text",x=-0.35,y=-0.50,parse=TRUE,size=5,label=paste('R:',NPS_Feces_Saliva_filter0.001.R1),family="serif",fontface="italic",colour="darkred")+annotate("text",x=-0.34,y=-0.47,parse=TRUE,size=5,label=paste('p:',NPS_Feces_Saliva_filter0.001.P1),family="serif",fontface="italic",colour="darkred")+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.position = "right",legend.title=element_text(size=5))
   dev.off()
  ##Figure1B
   pdf("HCNandNBT.PCoA-dist_20210331.pdf",width=3,height=3)
   HCNandNBT.PCoA<- ggplot(HCNandNBTdata,aes(x=x,y=y,color=HCNandNBTtype))+scale_color_manual(values=c('HCN'="#3F59A8",'NBT'="#EE3725"),labels=c('N_HC(N=44)','N_NPC(N=35)'))+
   geom_point(alpha=.7,size=2)+labs(x=paste("PCoA1(",format(100*HCNandNBTeig[1]/sum(HCNandNBTeig),digits=4),"%)",sep=""),y=paste("PCoA2(",format(100*HCNandNBTeig[2]/sum(HCNandNBTeig),digits=4),"%)",sep=""))+
   theme_classic()+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=6))+
   theme(legend.title = element_blank(), legend.position = c(0.8, 0.9)) +guides(fill=FALSE)+
   stat_ellipse(type="t",linetype=2)+annotate("text",x=0.40,y=-0.45,parse=TRUE,size=2,label=paste('R:',0.046),colour="#000000")+annotate("text",x=0.41,y=-0.50,parse=TRUE,size=2,label=paste('p:',HCNandNBTP1),colour="#000000")
   HCNandNBT.PCoA1<-ggMarginal(HCNandNBT.PCoA,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
   print(HCNandNBT.PCoA1)
   dev.off()
  ##Figure1C
   pdf("HCSandSBT.PCoA-dist_20210331.pdf",width=3,height=3)
   HCSandSBT.PCoA<- ggplot(HCSandSBTdata,aes(x=x,y=y,color=HCSandSBTtype))+scale_color_manual(values=c('HCS'="#3F59A8",'SBT'="#EE3725"),labels=c('S_HC(N=29)','S_NPC(N=37)'))+
   geom_point(alpha=.7,size=2)+labs(x=paste("PCoA1(",format(100*HCSandSBTeig[1]/sum(HCSandSBTeig),digits=4),"%)",sep=""),y=paste("PCoA2(",format(100*HCSandSBTeig[2]/sum(HCSandSBTeig),digits=4),"%)",sep=""))+
   theme_classic()+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=6))+
   theme(legend.position='none')+ ###去除整个legend
   stat_ellipse(type="t",linetype=2)+annotate("text",x=0.40,y=-0.27,parse=TRUE,size=2,label=paste('R:',0.107),colour="#000000")+annotate("text",x=0.41,y=-0.31,parse=TRUE,size=2,label=paste('p:',HCSandSBTP1),colour="#000000")
   HCSandSBT.PCoA1<-ggMarginal(HCSandSBT.PCoA,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
   print(HCSandSBT.PCoA1)
   dev.off()
 ###Figure1D
  pdf("HCFandFBT.PCoA-dist_20210331.pdf",width=3,height=3)
   HCFandFBT.PCoA<- ggplot(HCFandFBTdata,aes(x=x,y=y,color=HCFandFBTtype))+scale_color_manual(values=c('HCF'="#3F59A8",'FBT'="#EE3725"),labels=c('F_HC(N=15)','F_NPC(N=35)'))+
   geom_point(alpha=.7,size=2)+labs(x=paste("PCoA1(",format(100*HCFandFBTeig[1]/sum(HCFandFBTeig),digits=4),"%)",sep=""),y=paste("PCoA2(",format(100*HCFandFBTeig[2]/sum(HCFandFBTeig),digits=4),"%)",sep=""))+
   theme_classic()+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=6))+
   theme(legend.position='none')+ ###去除整个legend
   stat_ellipse(type="t",linetype=2)+annotate("text",x=0.40,y=-0.37,parse=TRUE,size=2,label=paste('R:',0.035),colour="#000000")+annotate("text",x=0.41,y=-0.41,parse=TRUE,size=2,label=paste('p:',HCFandFBTP1),colour="#000000")
   HCFandFBT.PCoA1<-ggMarginal(HCFandFBT.PCoA,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
   print(HCFandFBT.PCoA1)
   dev.off()

###Figure2 --20211230
 cd /public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2
 conda activate r_env
 R
 ##Figure2A
  pdf("Figure2A_Saliva_Genus.adonis_20211230.pdf",width=2.5,height=2.5)
      ggplot(Saliva_Genus_FeatureID0.001_PCoA_R2_2_phylum, aes(x = genus, y = singleR2, fill=Phylum)) +geom_bar(stat = "identity") + 
      #scale_fill_manual(values=c("Saliva"="#228B22","Saliva"="#FFA500","Feces"="#6495ED")) +
      labs(x=paste("genus"),y=paste("Adonis R2"))+theme_classic()+scale_fill_manual(values=c("Firmicutes"="#C25E74","Proteobacteria"="#00C7FE","Actinobacteria"="#F9F871","Bacteroidetes"="#618456","Fusobacteria"="#69A7FF"))+
      theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1,size=7,face="italic",colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=5,face="italic"),legend.title = element_text(color ="black", size = 7))+
      theme(legend.position=c(0.8,0.7))+theme(legend.key.size = unit(10, "pt"))+theme(legend.background = element_rect(fill = "white", colour = "black", size = 0.3 ) )      
      dev.off()
 ##Figure2B
  library(factoextra)
   p<-fviz_pca_biplot(HCSandSBTgenus.sum.res1_filter0.001.t.pca, label = c("g__Streptococcus","g__Prevotella","g__Neisseria"), habillage=as.character(HCSandSBTgenus.sum.res1_filter0.001.t1_group1$group),
                addEllipses=TRUE, ellipse.level=0.95,win.asp=0.8,geom.ind="point",fill.ind=HCSandSBTgenus.sum.res1_filter0.001.t1_group1$group,col.ind="black",palette=c('HCS'="#3F59A8",'SBT'="#EE3725"),
                ggtheme = theme_minimal())
   ggsave("Figure2B_Saliva_PCA_20211230.pdf",p,width=2.5,height=2.5)
 ##Figure2C
  pdf("Figure2C_HCandNPC-BT-Streptococcus_Genus_20210417.pdf",width=2.5,height=2.5)
  ggplot(HCSandSBTgenus.sum.res1_filter0.001_2t2_melt_group_g__Streptococcus,aes(Treatmenttype,value*100,fill=group))+geom_boxplot(aes(fill=group),outlier.colour=NA)+
  scale_fill_manual(values=c('HCS'="#3F59A8",'SBT'="#EE3725"),labels=c('HCS(N=29)','SBT(N=37)'))+
  theme_classic()+labs(x=paste(""),y=paste("relative abundance (%)\n Streptococcus"),title="")+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))+geom_jitter(position=position_jitter(0.08),size=0.5,alpha=0.8,colour="#000000")+
  theme(legend.position='right')+
  geom_signif(comparisons=list(c('HCS','SBT')),test=wilcox.test,step_increase=0.1,map_signif_level=T, text=(size=2),colour="#000000")
  dev.off()
 ##Figure2D
  pdf("Figure2C_Saliva_HCS_Streptococcus_20211001.pdf",height=4,width=2.5)
  ggplot(data=NPC288metagenome_metaphlan_species_Saliva3t_Streptococcus_HC1_melt,aes(x=variable,y=value1))+geom_boxplot(aes(fill=color))+theme_classic()+labs(x=paste("Species"),y=paste("log10(relative abundance)"),title="HC")+scale_fill_manual(values=c("#3F59A8"))+
   theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,face="italic",colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))+
   coord_flip()+theme(legend.position="none")
  dev.off()

 ##Figure2E
   color<-c("#F9F871","#FFC45F","#00B579","#FF9174","#71C25E","#00C7FE","#F3EED9","#69A7FF","#DBA053","#A06C20")
  Saliva_Streptococcus_Species1.melt_lefse<-filter(Saliva_Streptococcus_Species1.melt,Species=="S.anginosus"|Species=="S.cristatus"|Species=="S.gordonii"|Species=="S.infantis"|Species=="S.parasanguinis"|Species=="S.peroris"|Species=="S.pseudopneumoniae"|Species=="S.sanguinis"|Species=="S.sp_SK140"|Species=="S.vestibularis")
  Saliva_Streptococcus_Species1_lefse<-filter(Saliva_Streptococcus_Species1,Species=="S.anginosus"|Species=="S.cristatus"|Species=="S.gordonii"|Species=="S.infantis"|Species=="S.parasanguinis"|Species=="S.peroris"|Species=="S.pseudopneumoniae"|Species=="S.sanguinis"|Species=="S.sp_SK140"|Species=="S.vestibularis")
  colnames(Saliva_Streptococcus_Species1_lefse)<-c("Species","HC","NPC")
  Saliva_Streptococcus_Species1.melt_lefse$variable<-gsub("Saliva_","",Saliva_Streptococcus_Species1.melt_lefse$variable)

  p<-ggplot(Saliva_Streptococcus_Species1.melt_lefse, aes(x=variable, y=value, fill=Species)) + 
  geom_bar(stat = "identity", width=0.5, col='black') + theme_classic()+scale_fill_manual(values=color)+theme(legend.text=element_text(size = 7,face="italic",colour="#000000"))+labs(x=paste(""),y=paste("relative abundance(%)"))+
  geom_segment(data=Saliva_Streptococcus_Species1_lefse %>% arrange(by=desc(Species)) %>% mutate(HC=cumsum(HC)) %>% mutate(NPC=cumsum(NPC)), aes(x=1.25, xend=1.75, y=HC, yend=NPC))  
  ggsave("Figure2D_Saliva_Streptococcus_Species_20211001.pdf",p,width=3.5, height =4)

 ##Figure2F
  cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome/metagenome_reanalysis
  conda activate MicrobiotaProcess
  conda activate r_env
  R
  ###HCS
  p_HCS <- ggtree(HCS_network_metaphlan_species_trda, layout="inward_circular", size=0.2, xlim=c(18,NA))
  p_HCS <- p_HCS %<+%  HCS_network_metaphlan_species_melt_node3

  p_HCS1 <- p_HCS +geom_tippoint(mapping=aes(color=Taxa,shape=Level),size=1,alpha=1) +scale_color_manual(values=c('Firmicutes'='#C0BCDB','Proteobacteria'='#E35444','Actinobacteria'='#328774','Fusobacteria'='#FDDB47','Bacteroidetes'='#0000FF'),
        guide=guide_legend(keywidth=0.5,keyheight=0.5,order=2,override.aes=list(shape=c("Firmicutes" =18,"Proteobacteria" =20,"Actinobacteria"=20,"Fusobacteria" =20,"Bacteroidetes" =20),size=3),na.translate=TRUE)) +
        scale_shape_manual(values=c("Phylum"=20), guide="none" )  ###根据门属种的聚类进行画聚类树

  p_HCS2 <- p_HCS1 +new_scale_color() +geom_taxalink(data=HCS_weight_edge2,mapping=aes(taxa1=Species1,taxa2=Species2,color=Interaction),alpha=1,offset=0.1,size=0.15,ncp=10,hratio=1,arrow=grid::arrow(length = unit(0.005, "npc"))) +
        scale_colour_manual(values=c("Positive"="#DA4B94", "Negative"="#009999"),guide=guide_legend(keywidth=0.8, keyheight=0.5,order=1, override.aes=list(alpha=1, size=0.5)))

  #p3 <- p2 +geom_fruit(data=BGCsda,geom=geom_star,mapping=aes(y=Strain,x=BGCs,size=Count,fill=BGCs),starshape = 13,starstroke = 0,offset=-0.9,pwidth=1.4,grid.params=list(linetype=3)) +
        #scale_size_continuous(range=c(0, 2),limits=c(sort(BGCsda$Count)[2], max(BGCsda$Count)),breaks=c(1, 2, 3),name=bquote(paste(Log[2],"(",.("Count+1"), ")")),guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=4,override.aes = list(starstroke=0.3))) +
        #scale_fill_manual(values=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3","#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=3))
 
  #p4 <- p2 +geom_tiplab(mapping=aes(label=names),align=TRUE,size=5,linetype=NA,offset=5)  ###offset调整两圈之间的距离

  p_HCS5 <- p_HCS2 + new_scale_fill() +geom_fruit(data=HCS_edgs_weight_melt,geom=geom_bar,mapping=aes(x=value,y=Species,fill=variable),stat="identity",orientation="y",offset=0.01,pwidth=1.2,axis.params=list(axis = "x",text.angle = -45,hjust = 0,vjust = 0.5,nbreak =6,size=7)) +
        scale_fill_manual(name = "Number of interactions",values=c("positive"="#DA4B94", "negative"="#009999"),guide=guide_legend(keywidth=0.5,keyheight=0.5,order=3)) +
        theme(legend.background=element_rect(fill=NA),legend.title=element_text(size=7),legend.text=element_text(size=7),legend.spacing.y = unit(0.02, "cm"),legend.margin=margin(0.1, 0.3, 0.1,-0.3, unit="cm"),legend.box.margin=margin(0.1, 0.9, 0.1, -0.9, unit="cm"),plot.margin = unit(c(-1.2, -1.2, -1.2, 0.1),"cm"))
  ggsave("HCS_Species_network_20211022.pdf",p_HCS5,width=4,height=3)

 ##SBT
    pSBT <- ggtree(SBT_network_metaphlan_species_trda, layout="inward_circular", size=0.2, xlim=c(18,NA))
  pSBT <- pSBT %<+%  SBT_network_metaphlan_species_melt_node3

  pSBT1 <- pSBT +geom_tippoint(mapping=aes(color=Taxa,shape=Level),size=1,alpha=1) +scale_color_manual(values=c('Firmicutes'='#C0BCDB','Proteobacteria'='#E35444','Actinobacteria'='#328774','Fusobacteria'='#FDDB47','Bacteroidetes'='#0000FF',"Synergistetes"="#7FBC41","Spirochaetes"="#4D9221","Candidatus_Saccharibacteria"="#276419"),
        guide=guide_legend(keywidth=0.5,keyheight=0.5,order=2,override.aes=list(shape=c("Firmicutes" =18,"Proteobacteria" =20,"Actinobacteria"=20,"Fusobacteria" =20,"Bacteroidetes" =20,"Synergistetes"=20,"Spirochaetes"=20,"Candidatus_Saccharibacteria"=20),size=3),na.translate=TRUE)) +
        scale_shape_manual(values=c("Phylum"=20), guide="none" )  ###根据门属种的聚类进行画聚类树

  pSBT2 <- pSBT1 +new_scale_color() +geom_taxalink(data=SBT_weight_edge2,mapping=aes(taxa1=Species1,taxa2=Species2,color=Interaction),alpha=1,offset=0.1,size=0.15,ncp=10,hratio=1,arrow=grid::arrow(length = unit(0.005, "npc"))) +
        scale_colour_manual(values=c("Positive"="#DA4B94", "Negative"="#009999"),guide=guide_legend(keywidth=0.8, keyheight=0.5,order=1, override.aes=list(alpha=1, size=0.5)))

  #p3 <- p2 +geom_fruit(data=BGCsda,geom=geom_star,mapping=aes(y=Strain,x=BGCs,size=Count,fill=BGCs),starshape = 13,starstroke = 0,offset=-0.9,pwidth=1.4,grid.params=list(linetype=3)) +
        #scale_size_continuous(range=c(0, 2),limits=c(sort(BGCsda$Count)[2], max(BGCsda$Count)),breaks=c(1, 2, 3),name=bquote(paste(Log[2],"(",.("Count+1"), ")")),guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=4,override.aes = list(starstroke=0.3))) +
        #scale_fill_manual(values=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3","#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=3))
 
 # p4 <- p2 +geom_tiplab(mapping=aes(label=names),align=TRUE,size=5,linetype=NA,offset=5)  ###offset调整两圈之间的距离

  pSBT5 <- pSBT2 + new_scale_fill() +geom_fruit(data=SBT_edgs_weight_melt,geom=geom_bar,mapping=aes(x=value,y=Species,fill=variable),stat="identity",orientation="y",offset=0.01,pwidth=1.2,axis.params=list(axis = "x",text.angle = -45,hjust = 0,vjust = 0.5,nbreak =6)) +
        scale_fill_manual(name = "Number of interactions",values=c("positive"="#DA4B94", "negative"="#009999"),guide=guide_legend(keywidth=0.5,keyheight=0.5,order=3)) +
        theme(legend.background=element_rect(fill=NA),legend.title=element_text(size=7),legend.text=element_text(size=7),legend.spacing.y = unit(0.02, "cm"),legend.margin=margin(0.1, 0.1, 0.1,-0.1, unit="cm"),legend.box.margin=margin(0.1, 0.9, 0.1, -0.9, unit="cm"),plot.margin = unit(c(-1.2, -1.2, -1.2, 0.1),"cm"))
  ggsave("SBT_Species_network_20211022.pdf",pSBT5,width=4.5,height=3)

 ###Figure2G
  
###Figure3
  cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome/metagenome_reanalysis
  conda activate r_env
  R
 ##Figure3A
 
  
 ##Figure3B
   p<-ggplot(importance_Saliva_species_Streptococcus_lefse1, aes(x =species , y = MeanDecreaseAccuracy)) + #geom_hline(yintercept = -4, color = "black", size = 1) + # 添加y=0的辅助线
    geom_point(aes(color = group), size = 2) + scale_color_manual(values=c("significant"="#CD3333","non-significant"="#8B8378"))+        # 将点的size设置大一些比较好看
    geom_bar(aes(fill = group), stat = "identity", width = 0.05) + # 注意将width宽度设小
    theme_bw(base_family = "Times") + scale_fill_manual(values=c("significant"="#CD3333","non-significant"="#8B8378")) +
    theme(panel.grid.minor = element_blank(),
       panel.grid.major.x = element_blank(),      # 消除竖条的背景线
        axis.text.x = element_text(angle = 90),
        legend.position = "None",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"), # Mac 电脑上绘图展现中文需要此行命令
        plot.title = element_text(hjust = 0.5)) +  # 标题居中，若无标题可不加
        labs(x = "species", y = "MeanDecreaseAccuracy",
       colour = "", linetype = "", fill = "")+
       coord_flip()+
   theme(legend.key=element_blank(), 
   axis.text.x = element_text(colour = "black", size = 7), 
   axis.text.y = element_text(colour = "black", face = "italic", size = 7), 
   legend.text = element_text(size = 7, colour ="black"), 
   legend.title = element_text(size = 7, face = "bold"), 
   panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
   legend.position = "none") +ylim(-4,16)
  ggsave("Saliva_Streptococcus_species_importance_20211223.pdf",p,width=4,height=3)

 ###Figure3C
  pdf("Figure3C_HCS_SBT_Streptococcus_20211231.pdf",width=6,height=3)
  ggplot(HCSandSBT_streptococcus_species3_melt_group,aes(group,value,fill=group))+geom_boxplot(aes(fill=group),outlier.colour=NA)+
  scale_fill_manual(values=c('S_HC'="#3F59A8",'S_NPC'="#EE3725"),labels=c('HCS(N=29)','SBT(N=37)'))+
  theme_classic()+labs(x=paste(""),y=paste("Relative abundance(%)"),title="")+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))+geom_jitter(position=position_jitter(0.08),size=0.5,alpha=0.8,colour="#000000")+
  theme(legend.position='right')+
  facet_wrap(~variable,scales="free")+theme(strip.text.x=element_text(size=7))+geom_signif(comparisons=list(c('S_HC','S_NPC')),test=wilcox.test,step_increase=0.1,map_signif_level=T, text=(size=2),colour="#000000")
  dev.off()
  
 ##Figure3D
   pdf("Figure3D_Saliva.EBV_Streptococcus_ROC_20211024.pdf",width=4,height=4)
    par(lty=1)
    plot(EBV_copy20210414_Saliva1.fit.full.pre.data.roc,print.auc=FALSE,printz.thres=TRUE,print.thres.cex=0.9,print.auc.x=0.3,print.auc.y=0.4,auc.polygon=TRUE,auc.polygon.col="white",grid=c(0.5,0.2),grid.col=c("white","white"),max.auc.polygon=FALSE,smooth=T,main="",
    col= "black",legacy.axes=TRUE,cex.main=0.5,cex.axis=0.5,cex.lab=0.5,lwd=1)
    plot(NPC288metagenome_metaphlan_species_Saliva3t_Streptococcus_lefse3_EBV2.fit.full.pre.data.roc,add=T,col="red",print.auc=FALSE,print.auc.x=0.3,print.auc.y=0.3,smooth=T,lwd=1) 
    legend(0.80,0.25,bty="n",title="",legend=c("EBV(AUC:0.66)","Streptococcus_EBV(AUC:0.86)"),col=c("black","red"),lwd=1,cex=0.4)
    text(0.4, 0.3, labels=paste("P value =0.02"), adj=c(0, .5))
    dev.off()

 ##Figure3E
     mytheme <- theme_classic()+theme(axis.title=element_text(size=7,colour = 'black'), #坐标轴标题
                  axis.text.x = element_text(colour = "black", size = 7), axis.text.y = element_text(colour = "black", face = "italic", size = 7),  #坐标轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank(),
                 legend.text=element_text(size=7),#关闭图例边框
                 legend.position="none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major= element_blank())

   Saliva_Streptococcus_species_EBV_ROC1$ymin<-rep(0.6,nrow(Saliva_Streptococcus_species_EBV_ROC1))
   Saliva_Streptococcus_species_EBV_ROC1$names<-gsub("Streptococcus","S.(3 species)",Saliva_Streptococcus_species_EBV_ROC1$names)
   Saliva_Streptococcus_species_EBV_ROC1$names<-factor(Saliva_Streptococcus_species_EBV_ROC1$names,levels=rev(c("S.(3 species)_EBV","S.(3 species)","S.infantis_EBV","S.infantis","S.gordonii_EBV","S.gordonii","S.sanguinis_EBV","S.sanguinis","EBV")))

   p<-ggplot(Saliva_Streptococcus_species_EBV_ROC1,aes(names,ROC)) +
   geom_point(aes(fill=names,size=5),alpha=0.6,pch=21) +  #fill对应点的填充色，colour对应点的边框色
   geom_segment(aes(x=names,y=ymin,xend=names,yend=ROC))+
   scale_fill_manual(values=c("S.(3 species)_EBV"='#232F95', "S.(3 species)"='#8D002E',"S.infantis_EBV"='#232F95', "S.infantis"='#8D002E',"S.gordonii_EBV"='#232F95', "S.gordonii"='#8D002E',"S.sanguinis_EBV"='#232F95', "S.sanguinis"='#8D002E',"EBV"="#000000")) + #设定颜色的变化范围
   #scale_size_area(max_size = 8, breaks=c(2,4,6,8,10)) + #设定点的大小比例和图例上显示的间隔
   labs(y='ROC',x='')+ylim(0.6,0.93)+
   coord_flip()+mytheme
  ggsave("Figure3E_Saliva_Streptococcus_species_EBV_ROC_20211024.pdf",p,width=4,height=4)

###Figure4
  cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome
  conda activate r_env
  ##Figure 4A
   pdf("Figure4A_HCSandSBT.PCoA__KO_dist_20211203.pdf",width=3,height=3)
   HCSandSBT_KO_PCoA<- ggplot(HCSandSBT_KO_data,aes(x=x,y=y,color=HCSandSBT_KO_type))+scale_color_manual(values=c('HCS'="#3F59A8",'SBT'="#EE3725"),labels=c('HC(N=27)','NPC(N=35)'))+
   geom_point(alpha=.9,size=2)+labs(x=paste("PCoA1(",format(100*HCSandSBT_KO_eig[1]/sum(HCSandSBT_KO_eig),digits=4),"%)",sep=""),y=paste("PCoA2(",format(100*HCSandSBT_KO_eig[2]/sum(HCSandSBT_KO_eig),digits=4),"%)",sep=""))+
   theme_classic()+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=6))+
   theme(legend.title = element_blank(), legend.position = c(0.85, 0.20)) +
   #guides(fill=FALSE)+
   #theme(legend.position=c(0.85,0.20))+
   #scale_color_manual(values=c('HCF'="#ff0000",'FBT'="#0000FF"),labels=c('HCF(N=15)','FBT(N=35)'))
   stat_ellipse(type="t",linetype=2)+annotate("text",x=0.15,y=-0.20,parse=TRUE,size=2,label=paste('R:',0.252),colour="#000000")+annotate("text",x=0.27,y=-0.20,parse=TRUE,size=2,label=paste('p:',HCSandSBT_KO_P1),colour="#000000")
   #HCSandSBT_KO_PCoA1<-ggMarginal(HCSandSBT_KO_PCoA,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
   print(HCSandSBT_KO_PCoA)
   dev.off()
   
  ##Figure 4B
   cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome
   conda activate r_env
  
  pdf("Figure4B__S_HCandNPC_enrichKO_20220101.pdf",width=6,height=6)
  ggplot(S_HCandNPC_enrichKO_df3,aes(x=group,y=Description))+
  geom_point(aes(size=`GeneCount`,
                 color=`Pathway1`))+
  theme_bw()+
  scale_color_manual(values=c("Metabolism"="#69A7FF","Genetic Information Processing"="#FF6997","Brite Hierarchies"="#618456","Environmental Information Processing"="#2F4858","Cellular Processes"="#CDAD00","Human Diseases"="#C25E74"))+
  labs(x=paste("Group"),y=paste("KEGG Pathway"))+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=45,hjust = 0.5,vjust=0.6,size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7),legend.title = element_text(color ="black", size = 7))+
  coord_flip()+theme(legend.position='bottom')
  dev.off()
 
 ###Figure4D
  cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome
  conda activate r_env
  pdf("Saliva_Species_Pathway3_spearman_20220103.pdf",width = 5,height = 5)  ##pdf(file = "heatmap.pdf",width = 5,height = 5,family="Times")
  col3 <- colorRampPalette(c("white","#00008B")) 
  corrplot(Saliva_Species_Pathway3_cor_dcast4, bg="white",title=" ",type= "full", tl.srt = 60,is.corr = FALSE,col=col3(10), outline = FALSE,tl.pos = NULL,tl.cex = 0.6, tl.col = "black",mar = c(1, 1, 1, 1), method="circle",addgrid.col = NULL)
  
   dev.off()
 ###Figure4E
    p<-ggplot(S_Pathway3_metaphlan_species_merge_Streptococcus_mitis_Two_component_system,aes(x=s__Streptococcus_mitis_oralis_pneumoniae,y=Two_component_system*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=S_Pathway3_metaphlan_species_merge_Streptococcus_mitis_Two_component_system,method="spearman",geom="text")+
    theme_classic()+labs(x=paste("Streptococcus mitis"),y=paste("relative abundance(%)\n Two_component_system"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure4D_Streptococcus mitis_Two_component_system_20211205.pdf",p,width=2,height=2)

###Figure5
 ###Figure 5A
   cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome/metagenome_reanalysis
   pdf("Figure5A_Saliva_EBV_meta_Streptococcus_species.correlation_20211004.pdf",width=4,height=5)
  ggplot(data=Saliva_EBV_meta_bacteria_na.omit2_z_score_cor_select2.melt2,aes(x=variable,y=species))+
  geom_point(aes(size=abs(value_1),
                 color=factor(value_1)),
             shape=15)+
  scale_color_manual(values = c(rep("#009ccc",4),rep("#fe0000",4)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="grey"),
        axis.ticks = element_blank())+
  geom_segment(data=Saliva_EBV_meta_bacteria_na.omit2_z_score_cor_select2.melt3,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey")+
  geom_segment(data=Saliva_EBV_meta_bacteria_na.omit2_z_score_cor_select2.melt4,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey")+
  scale_size_continuous(range = c(1,6))+
  scale_y_discrete(position = "right")+
  labs(x=NULL,y=NULL)+
  theme(axis.text
 
 ###Figure 5B
  pdf("Figure5B_Saliva_Streptococcus_Clinical_summary_20211004.pdf",width=4,height=4)
   ggplot(Saliva_Streptococcus_clinical.melt, aes(x=variable,y=value))+geom_bar(aes(fill=group),stat="identity")+scale_fill_lancet()+
   coord_polar(theta = "x")+ylim(0,20)+theme_light()+labs(x=paste(),y=paste("Number of associations"))
   dev.off()

   p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=NLR,y=s__Streptococcus_gordonii*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("NLR"),y=paste("relative abundance(%)\n S.gordonii"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_NLR_S.gordonii_20211004.pdf",p,width=2.5,height=2.5)

  p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=NLR,y=s__Streptococcus_cristatus*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("NLR"),y=paste("relative abundance(%)\n S.cristatus"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_NLR_S.cristatus_20211004.pdf",p,width=2.5,height=2.5) ###有意义


   p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=NLR,y=s__Streptococcus_mitis_oralis_pneumoniae*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+#stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+
    theme_classic()+ labs(x=paste("NLR"),y=paste("relative abundance(%)\n S.mitis_oralis_pneumoniae"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_NLR_S.mitis_oralis_pneumoniae_20211004.pdf",p,width=2.5,height=2.5)

   p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=NLR,y=s__Streptococcus_mutans*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("NLR"),y=paste("relative abundance(%)\n S.mutans"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_NLR_S.mutans_20211004.pdf",p,width=2.5,height=2.5)


   p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=NLR,y=s__Streptococcus_infantis*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("NLR"),y=paste("relative abundance(%)\n S.infantis"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_NLR_S.infantis_20211004.pdf",p,width=2.5,height=2.5)

   p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=NLR,y=s__Streptococcus_sanguinis*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("NLR"),y=paste("relative abundance(%)\n S.sanguinis"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_NLR_S.sanguinis_20211004.pdf",p,width=2.5,height=2.5)



    p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=PLR,y=s__Streptococcus_infantis*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("PLR"),y=paste("relative abundance(%)\n S.infantis"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5D_Saliva_PLR_S.infantis_20211004.pdf",p,width=2.5,height=2.5)

    p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=PLR,y=s__Streptococcus_gordonii*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("PLR"),y=paste("relative abundance(%)\n S.gordonii"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_PLR_S.gordonii_20211004.pdf",p,width=2.5,height=2.5)


     p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=PLR,y=s__Streptococcus_sp_oral_taxon_071*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("PLR"),y=paste("relative abundance(%)\n S.sp_oral_taxon_071"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_PLR_S.sp_oral_taxon_071_20211004.pdf",p,width=2.5,height=2.5)

   p<-ggplot(Saliva_EBV_meta_Streptococcus_species_na.omit1,aes(x=LMR,y=s__Streptococcus_sanguinis*100))+geom_point(color="red",size=1)+
    stat_smooth(method="lm",se=TRUE)+stat_cor(data=Saliva_EBV_meta_Streptococcus_species_na.omit1,method="spearman",geom="text")+theme_classic()+
    labs(x=paste("LMR"),y=paste("relative abundance(%)\n S.sanguinis"))+
    theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
   ggsave("Figure5C_Saliva_LMR_S.sanguinis_20211004.pdf",p,width=2.5,height=2.5)


###control--20210936
 cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome/metagenome_reanalysis
 conda activate r_env

 NPC314otu<-read.table("/public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2/NPC311otu.txt",header=T,check.names=FALSE)
 NPC314tax<-read.table("/public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2/NPC311taxonomy.txt",header=T,fill=T,check.names=FALSE)
 NPC314otu.tax<-merge(NPC314tax,NPC314otu,by="FeatureID")
 Control_otu.tax<-subset(NPC314otu.tax,select=c(FeatureID,Kingdom,Phylum,Class,Order,Family,Genus,Species,ControlM1,ControlN1,ControlP1,NPC002CN))
 Controlgenus<-subset(Control_otu.tax,select=-c(FeatureID,Kingdom,Phylum,Class,Order,Family,Species)) #保留genus水平
 dim(Controlgenus)
 library(dplyr)
 Controlgenus1 <-aggregate(Controlgenus[,2:ncol(Controlgenus)],list(Controlgenus[,1]),sum)  #合并相同的genus
 dim(Controlgenus1)
 Controlgenus2 <- data.frame(Controlgenus1[-1,],check.names=FALSE) #去掉空白行
 Controlgenus2 <- data.frame(Controlgenus2[-1,],check.names=FALSE) #去掉g__行
 dim(Controlgenus2) #663个属
 names(Controlgenus2)[names(Controlgenus2)=='Group.1']<-'Genus'
 Controlgenus3<- mutate(Controlgenus2,sum=apply(Controlgenus2[,2:ncol(Controlgenus2)],1,sum))
 Controlgenus3 <- Controlgenus3[order(Controlgenus3$sum,decreasing=T),]
 Controlgenus3<-filter(Controlgenus3, sum >0) #去掉read数总和小于100的属
 dim(Controlgenus3) #136个属
 Controlgenus4 <- subset(Controlgenus3,select=-c(sum))
 colSums(Controlgenus4[,(2:5)])
 rownames(Controlgenus4)<-Controlgenus4$Genus
 Controlgenus5<-subset(Controlgenus4,select=-c(Genus))
 Controlgenus5.res<-as.data.frame(t(t(Controlgenus5)/colSums(Controlgenus5)))
 Controlgenus5.res2<-subset(Controlgenus5.res,select=c(ControlM1,ControlP1))
 Controlgenus5.res3<-tibble::rownames_to_column(Controlgenus5.res2,"Genus")
 Controlgenus5.res3_Streptococcus<-filter(Controlgenus5.res3,Genus=="g__Streptococcus")
 colnames(Controlgenus5.res3_Streptococcus)<-c("Genus","Control1","Control2")
 Controlgenus5.res3_Streptococcus.melt<-melt(Controlgenus5.res3_Streptococcus,id.vars="Genus")
 Controlgenus5.res3_Streptococcus.melt$value<-as.numeric(as.character(Controlgenus5.res3_Streptococcus.melt$value))

  pdf("Figure_control_HCandNPC-BT_Streptococcus_20210926.pdf",width=3,height=4)
  ggplot(Controlgenus5.res3_Streptococcus.melt,aes(variable,value*100,fill=variable))+geom_bar(stat = "identity",position=position_dodge(0.9))+
  scale_fill_manual(values=c('Control1'="#3F59A8",'Control2'="#EE3725"),labels=c('Control1','Control2'))+
  theme_classic()+labs(x=paste(""),y=paste("relative abundance (%)"),title="")+theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1,size=0,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))
  dev.off()
###Figure S1
  cd /public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2
 
 conda activate r_env
 R
 sampleType.value$value<-as.numeric(as.character(sampleType.value$value))

 pdf("Figure_S1C_sampleType.adonis_20220106.pdf",width=4,height=3)
 ggplot(sampleType.value, aes(x = sampleType, y = value, fill=sampleType)) +geom_bar(stat = "identity") + scale_fill_manual(values=c('Feces'="#6495ED","Saliva"="#FFA500","NPS"="#228B22")) +labs(x=paste("sampleType"),y=paste("Adonis R2"))+theme_classic()+theme(axis.text.x=element_text(size=10,colour="#000000"),axis.text.y=element_text(size=10,colour="#000000"),axis.title.x=element_text(size=12,face="bold",colour="#000000"),axis.title.y=element_text(size=12,face="bold",colour="#000000"),legend.text=element_text(size=10))
 dev.off()

###Figure S2

 pdf("Figure_S2A_HCSandSBT.PCoA-dist_Streptococcus_20220106.pdf",width=3,height=3)
      HCSandSBT_Streptococcus.PCoA<- ggplot(HCSandSBTdata_Streptococcus,aes(x=x,y=y,color=HCSandSBTtype))+scale_color_manual(values=c('HCS'="#3F59A8",'SBT'="#EE3725"),labels=c('S_HC(N=29)','S_NPC(N=37)'))+
      geom_point(alpha=.9,size=2)+labs(x=paste("PCoA1(",format(100*HCSandSBTeig_Streptococcus[1]/sum(HCSandSBTeig_Streptococcus),digits=4),"%)",sep=""),y=paste("PCoA2(",format(100*HCSandSBTeig_Streptococcus[2]/sum(HCSandSBTeig_Streptococcus),digits=4),"%)",sep=""))+
      theme_classic()+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=6))+
      theme(legend.position='none')+ ###去除整个legend
      stat_ellipse(type="t",linetype=2)+annotate("text",x=0.35,y=-0.27,parse=TRUE,size=2,label=paste('R:',0.088),colour="#000000")+annotate("text",x=0.36,y=-0.31,parse=TRUE,size=2,label=paste('p:',HCSandSBT_Streptococcus_P1),colour="#000000")
      HCSandSBT_Streptococcus.PCoA1<-ggMarginal(HCSandSBT_Streptococcus.PCoA,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
      print(HCSandSBT_Streptococcus.PCoA1)
 dev.off()

###Figure S3
  cd /public4/lym/HCandNPC/NPC_metagenome/HC_NPC_pre_metagenome/metagenome_reanalysis
  conda activate ROC
  R
 ###用Streptococcus_infantis为点区分HCS and SBT----20210924
  NPC288metagenome_metaphlan_species_Saliva3t_Streptococcus_infantis<-subset(NPC288metagenome_metaphlan_species_Saliva3t_Streptococcus,select=c(s__Streptococcus_infantis))
  NPC288metagenome_metaphlan_species_Saliva3t_Streptococcus_infantis$group<-c(rep(c("S_HC"),26),rep(c("S_NPC"),35))
  Saliva_species_Streptococcus_infantis<-NPC288metagenome_metaphlan_species_Saliva3t_Streptococcus_infantis
  Saliva_species_Streptococcus_infantis$group<-factor(Saliva_species_Streptococcus_infantis$group,levels=c("S_NPC","S_HC"))
  
  Saliva_species_Streptococcus_infantis1<-tibble::rownames_to_column(Saliva_species_Streptococcus_infantis,"sampleID")
  HCSandSBTdata1<-read.table("/public/home/lym/Project/NPC/clinicalsample/MGGZ20190822A1463B-4N_A20191129-3834/qiime2-NPC-clinicalsamples/NPC314-reananlysis/reanalysis2/HCSandSBTdata1.txt",header=T)
  HCSandSBTdata1_Streptococcus_infantis1<-merge(HCSandSBTdata1,Saliva_species_Streptococcus_infantis1,by="sampleID")

  pdf("HCSandSBT_Genus0.001.PCoA_Streptococcus_20210924.pdf",width=6,height=4)
   ggplot(HCSandSBTdata1_Streptococcus_infantis1,aes(x=x,y=y))+geom_point(aes(fill=s__Streptococcus_infantis,shape=HCSandSBTtype),colour="black",alpha=.8,size=2)+scale_fill_gradient2(low = '#ff0000', high = '#0000FF')+scale_shape_manual(values=c('HCS'=21,'SBT'=24))+
   labs(x=paste("PCoA1(32.27%)",sep=""),y=paste("PCoA2(16.62%)",sep=""),title="HCSandSBTPCoA")+theme_classic()+annotate("text",x=0.35,y=-0.35,parse=TRUE,size=4,label=paste('R:0.12'),family="serif",fontface="italic",colour="darkred")+annotate("text",x=0.34,y=-0.39,parse=TRUE,size=4,label=paste('p:0.001'),family="serif",fontface="italic",colour="darkred")
   dev.off()
###Figure S4
 ####处理taxonomy.tsv
 cd /public4/lym/SRA_download/NPC_oral_microbiome_msystem/qiime2/NPC_Oral_microbiome_taxonomy
 cd /public4/lym/SRA_download/NPC_oral_microbiome_msystem/qiime2/R
 NPC_Oral_microbiome_msystem_taxonomy<-read.delim("/public4/lym/SRA_download/NPC_oral_microbiome_msystem/qiime2/NPC_Oral_microbiome_taxonomy/NPC_Oral_microbiome_msystem_taxonomy4.txt")
 NPC_Oral_microbiome_otu_table<-read.delim("/public4/lym/SRA_download/NPC_oral_microbiome_msystem/qiime2/NPC_Oral_microbiome_taxonomy/NPC_Oral_microbiome_otu-table1.tsv")

 PRJEB37445_metadata_20220326<-read.table("/public4/lym/SRA_download/NPC_oral_microbiome_msystem/qiime2/NPC_Oral_microbiome_taxonomy/PRJEB37445_metadata_20220326.txt",header=T)

 NPC_Oral_tax<-NPC_Oral_microbiome_msystem_taxonomy
 NPC_Oral_otu<-NPC_Oral_microbiome_otu_table
 NPC_Oral_meta<-PRJEB37445_metadata_20220326

 ###合并属
 NPC_Oral_otu.tax <-merge(NPC_Oral_tax,NPC_Oral_otu)
 head(NPC_Oral_otu.tax)
 dim(NPC_Oral_otu.tax)
 NPC_Oral_genus<-subset(NPC_Oral_otu.tax,select=-c(Feature_ID,Kingdom,Phylum,Class,Order,Family,Species)) #保留genus水平
 dim(NPC_Oral_genus)
 library(dplyr)
 NPC_Oral_genus1 <-aggregate(NPC_Oral_genus[,2:ncol(NPC_Oral_genus)],list(NPC_Oral_genus[,1]),sum)  #合并相同的genus
 dim(NPC_Oral_genus1)
 NPC_Oral_genus2 <- data.frame(NPC_Oral_genus1[-1,],check.names=FALSE) #去掉空白行
 ##NPC_Oral_genus2 <- data.frame(NPC_Oral_genus2[-1,],check.names=FALSE) #去掉g__行
 dim(NPC_Oral_genus2) #254个属
 names(NPC_Oral_genus2)[names(NPC_Oral_genus2)=='Group.1']<-'Genus'
 NPC_Oral_genus3<- mutate(NPC_Oral_genus2,sum=apply(NPC_Oral_genus2 [,2:ncol(NPC_Oral_genus)],1,sum))
 NPC_Oral_genus3 <- NPC_Oral_genus3 [order(NPC_Oral_genus3$sum,decreasing=T),]
 NPC_Oral_genus3<-filter(NPC_Oral_genus3, sum >0) #去掉read数总和小于100的属
 dim(NPC_Oral_genus3) #254个属
 NPC_Oral_genus4 <- subset(NPC_Oral_genus3,select=-c(sum))

 NPC_Oral_genus4_filter<-NPC_Oral_genus4[rowSums(NPC_Oral_genus4 >= 1)>=ncol(NPC_Oral_genus4)/10, ] ###至少在1/10的样本中的reads数大于1
 rownames(NPC_Oral_genus4_filter)<-NPC_Oral_genus4_filter$Genus
 NPC_Oral_genus4_filter1<-subset(NPC_Oral_genus4_filter,select=-c(Genus))
 head(NPC_Oral_genus4_filter1)
 NPC_Oral_genus4_filter1.res<-as.data.frame(t(t(NPC_Oral_genus4_filter1)/colSums(NPC_Oral_genus4_filter1)))
 NPC_Oral_genus4_filter1.res.t<-as.data.frame(t(NPC_Oral_genus4_filter1.res))
 NPC_Oral_genus4_filter1.res.t1<-tibble::rownames_to_column(NPC_Oral_genus4_filter1.res.t,"sampleID")

 NPC_Oral_genus4_filter1.res.t1_group<-merge(NPC_Oral_meta,NPC_Oral_genus4_filter1.res.t1)
 NPC_Oral_genus4_filter1.res.t1_group$group<-factor(NPC_Oral_genus4_filter1.res.t1_group$group,levels=c("control","case"))

 pdf("NPC_Oral-Streptococcus_Genus_20210326_2.pdf",width=2.5,height=2.5)
  ggplot(NPC_Oral_genus4_filter1.res.t1_group,aes(group,g_Streptococcus*100,fill=group))+geom_boxplot(aes(fill=group),outlier.colour=NA)+
  scale_fill_manual(values=c('control'="#3F59A8",'case'="#EE3725"),labels=c('HC(N=496)','NPC(N=499)'))+
  theme_classic()+labs(x=paste(""),y=paste("relative abundance (%)\n Streptococcus"),title="")+theme(axis.text.x=element_text(size=7,colour="#000000"),axis.text.y=element_text(size=7,colour="#000000"),axis.title.x=element_text(size=7,colour="#000000"),axis.title.y=element_text(size=7,colour="#000000"),legend.text=element_text(size=7))+geom_jitter(position=position_jitter(0.08),size=0.5,alpha=0.8,colour="#000000")+
  theme(legend.position='right')+
  geom_signif(comparisons=list(c('control','case')),test=wilcox.test,step_increase=0.1,map_signif_level=T, text=(size=2),colour="#000000")
  dev.off()

 set.seed(123)
 NPC_Oral_genus4_filter1.res.t1_group1<-subset(NPC_Oral_genus4_filter1.res.t1_group,select=-c(sampleID,run_accession,age,gender))

 #为了方便后续评估随机森林模型的性能
 #将总数据集分为训练集（占 70%）和测试集（占 30%）
 set.seed(123)
 NPC_Oral_genus4_filter1.res.t1_group1_train <- sample(nrow(NPC_Oral_genus4_filter1.res.t1_group1), nrow(NPC_Oral_genus4_filter1.res.t1_group1)*0.7)
 NPC_Oral_genus4_filter1.res.t1_group1_train1 <- NPC_Oral_genus4_filter1.res.t1_group1[NPC_Oral_genus4_filter1.res.t1_group1_train, ]
 NPC_Oral_genus4_filter1.res.t1_group1_test <- NPC_Oral_genus4_filter1.res.t1_group1[-NPC_Oral_genus4_filter1.res.t1_group1_train, ]
 ##randomForest 包的随机森林
 library(randomForest)
 #随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
 set.seed(123)
 NPC_Oral_genus4_filter1.res.t1_group1_train.forest <- randomForest(group~., data = NPC_Oral_genus4_filter1.res.t1_group1_train1, importance = TRUE)
 NPC_Oral_genus4_filter1.res.t1_group1_train.forest

 #% Var explained体现了预测变量（用于回归的所有OTU）对响应变量（植物年龄）有关方差的整体解释率,可以理解为该回归的R2=0.8917，相当可观的一个数值，表明该植物根际细菌的群落结构随植物生长密切相关
 #使用训练集，查看预测精度
  NPC_Oral_genus4_filter1.res.t1_group1_predict <- predict(NPC_Oral_genus4_filter1.res.t1_group1_train.forest,  NPC_Oral_genus4_filter1.res.t1_group1_train1)
 plot(NPC_Oral_genus4_filter1.res.t1_group1_train1$group,  NPC_Oral_genus4_filter1.res.t1_group1_predict, main = 'Train', xlab = 'group', ylab = 'Predict')
 abline(1, 1)
 #使用测试集，评估预测性能
  NPC_Oral_genus4_filter1.res.t1_group1_predict <- predict( NPC_Oral_genus4_filter1.res.t1_group1_train.forest,  NPC_Oral_genus4_filter1.res.t1_group1_test)
 plot( NPC_Oral_genus4_filter1.res.t1_group1_test$group,  NPC_Oral_genus4_filter1.res.t1_group1_predict, main = 'Test',
    xlab = 'group', ylab = 'Predict')
 abline(1, 1)
 ###评估预测变量的重要性
 ##OTU 的重要性评估
 #查看表示每个预测变量（细菌 OTU）重要性的得分
 #summary(otu_train.forest)
 importance_NPC_Oral_genus4_filter1.res.t1_group1 <- NPC_Oral_genus4_filter1.res.t1_group1_train.forest$importance
 head(importance_NPC_Oral_genus4_filter1.res.t1_group1)
 #或者使用函数 importance()
 importance_NPC_Oral_genus4_filter1.res.t1_group1 <- data.frame(importance(NPC_Oral_genus4_filter1.res.t1_group1_train.forest), check.names = FALSE)
 head(importance_NPC_Oral_genus4_filter1.res.t1_group1)
 #作图展示 top30 重要的 OTUs
 varImpPlot(NPC_Oral_genus4_filter1.res.t1_group1_train.forest, n.var = min(30, nrow(NPC_Oral_genus4_filter1.res.t1_group1_train.forest$importance)),
    main = 'Top 30 - variable importance')
 ##“%IncMSE”即increase in mean squared error，通过对每一个预测变量随机赋值，如果该预测变量更为重要，那么其值被随机替换后模型预测的误差会增大。因此，该值越大表示该变量的重要性越大；
 ##“IncNodePurity”即increase in node purity，通过残差平方和来度量，代表了每个变量对分类树每个节点上观测值的异质性的影响，从而比较变量的重要性。该值越大表示该变量的重要性越大。
 ###对于“%IncMSE”或“IncNodePurity”，二选一作为判断预测变量重要性的指标。需注意的是，二者的排名存在一定的差异。
  #可以根据某种重要性的高低排个序，例如根据“IncNodePurity”指标
 importance_NPC_Oral_genus4_filter1.res.t1_group2 <- importance_NPC_Oral_genus4_filter1.res.t1_group1[order(importance_NPC_Oral_genus4_filter1.res.t1_group1$MeanDecreaseAccuracy, decreasing = TRUE), ]
 head(importance_NPC_Oral_genus4_filter1.res.t1_group2)
 #输出表格
 write.table(importance_NPC_Oral_genus4_filter1.res.t1_group1, 'importance_NPC_Oral_genus4_filter1.res.t1_group1.txt', sep = '\t', col.names = NA, quote = FALSE)

 ###画图展示：

     mytheme <- theme_classic()+theme(axis.title=element_text(size=7,colour = 'black'), #坐标轴标题
                  axis.text.x = element_text(colour = "black", size = 7), axis.text.y = element_text(colour = "black", face = "italic", size = 7),  #坐标轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank(),
                 legend.text=element_text(size=7),#关闭图例边框
                 legend.position="none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major= element_blank())
 
 importance_NPC_Oral_genus4_filter1.res.t1_group1_top10<-importance_NPC_Oral_genus4_filter1.res.t1_group1[(1:10),]
 library(ggplot2)
 importance_NPC_Oral_genus4_filter1.res.t1_group1_top10_1<-tibble::rownames_to_column(importance_NPC_Oral_genus4_filter1.res.t1_group1_top10,"Genus")
 importance_NPC_Oral_genus4_filter1.res.t1_group1_top10_1$ymin<-rep(0.6,nrow(importance_NPC_Oral_genus4_filter1.res.t1_group1_top10_1))
 importance_NPC_Oral_genus4_filter1.res.t1_group1_top10_1$Genus<-factor(importance_NPC_Oral_genus4_filter1.res.t1_group1_top10_1$Genus,levels=rev(c("g_Streptococcus","g_Porphyromonas","g_Rothia","g_Fusobacterium","g_Neisseria","g_Veillonella","g_Haemophilus","g_Prevotella_7","g_Actinomyces","g_Alloprevotella")))


   p<-ggplot(importance_NPC_Oral_genus4_filter1.res.t1_group1_top10_1,aes(Genus,MeanDecreaseAccuracy)) +
   geom_point(aes(fill=Genus,size=5),alpha=0.6,pch=21) +  #fill对应点的填充色，colour对应点的边框色
   geom_segment(aes(x=Genus,y=ymin,xend=Genus,yend=MeanDecreaseAccuracy))+
   scale_fill_manual(values=c("g_Streptococcus"='#8D002E', "g_Porphyromonas"='#232F95',"g_Rothia"='#232F95', "g_Fusobacterium"='#232F95',"g_Neisseria"='#232F95', "g_Veillonella"='#232F95',"g_Haemophilus"='#232F95', "g_Prevotella_7"='#232F95',"g_Actinomyces"="#232F95","g_Alloprevotella"="#232F95")) + #设定颜色的变化范围
   #scale_size_area(max_size = 8, breaks=c(2,4,6,8,10)) + #设定点的大小比例和图例上显示的间隔
   labs(y='MeanDecreaseAccuracy',x='')+ylim(0.6,5.5)+
   coord_flip()+mytheme
  ggsave("NPC_Oral_Streptococcus_Genus_MeanDecreaseAccuracy_20220326_2.pdf",p,width=3,height=3)


 ###进行拟合回归
 ###将EBV，因此现拟合logistic回归模型
 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus<-subset(NPC_Oral_genus4_filter1.res.t1_group1,select=c(group,g_Streptococcus))

 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_HC<-filter(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus,group=="0")
 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_NPC<-filter(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus,group=="1")
 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_HC_mean<-colSums(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_HC)/nrow(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_HC)
 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_NPC_mean<-colSums(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_NPC)/nrow(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus_NPC)

 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group<-gsub("control","0",NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group)
 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group<-gsub("case","1",NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group)
 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group<-as.numeric(as.character(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group))

 NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full<-glm(group~.,data=NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus,family=binomial())
   summary(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full)
  NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre<-predict(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full,type='response')
  NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre.data<-data.frame(pro=NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre,obs= NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group)
  NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre.data.roc<-roc(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus$group, NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre)
   NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre.data.roc
   pROC::auc(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre.data.roc);pROC::ci(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre.data.roc) ##0.6645

 pdf("NPC_Oral_Streptococcus_ROC_20220326.pdf",width=3,height=3)
    par(lty=1)
    plot(NPC_Oral_genus4_filter1.res.t1_group1_Streptococcus.fit.full.pre.data.roc,print.auc=FALSE,printz.thres=TRUE,print.thres.cex=0.9,print.auc.x=0.3,print.auc.y=0.4,auc.polygon=TRUE,auc.polygon.col="white",grid=c(0.5,0.2),grid.col=c("white","white"),max.auc.polygon=FALSE,smooth=T,main="",
    col= "black",legacy.axes=TRUE,cex.main=0.5,cex.axis=0.5,cex.lab=0.5,lwd=1)
    legend(0.80,0.25,bty="n",title="",legend=c("Streptococcus(AUC:0.61)"),col=c("red"),lwd=1,cex=0.4)
    dev.off()


 https://github.com/lym20220107/lym20220107.git


