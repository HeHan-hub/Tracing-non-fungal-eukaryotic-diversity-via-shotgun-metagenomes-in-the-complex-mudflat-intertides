library(vegan)
library(ape)
library(iCAMP)
library(parallel)
library(bigmemory)
library(NST) 
library(permute)
library(DirichletReg)
library(linkET)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggcor)
library(dplyr)
setwd("E:/time-meta")

#Figgure 5a NST ratio
comm <- t(read.delim("eukspecies-time-normalised.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
group=read.delim("group.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

samp.ck=NST::match.name(rn.list=list(comm=comm,group=group))
comm=samp.ck$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
group=samp.ck$group

groupi=group[,1,drop=FALSE]
meta.groupi=groupi
dist.method="bray"

tnst=tNST(comm=comm, group=groupi, meta.group=meta.groupi, meta.com=NULL,
          dist.method=dist.method, abundance.weighted=TRUE, rand=1000,
          output.rand=TRUE, nworker=28, LB=FALSE, null.model="PF",dirichlet=F,
          between.group=T, SES=TRUE, RC=TRUE)

write.table(tnst$index.grp, file = "euk.tNST.summary.csv",sep =",", quote =FALSE)
write.table(tnst$index.pair.grp, "euk.tNST.pairwise.csv",quote = FALSE,sep = ",")
write.table(tnst$index.pair, "euk.tNST.pair.csv",quote = FALSE,sep = ",")

data<-read.csv("euk.tNST.csv")
data$Group <- factor(data$Group,levels=c("Aug","Oct","Dec","Feb","April","June"))
ggplot(data,aes(x=Group,y=NST))+
  stat_summary(aes(fill="Group"), alpha=0.65,fun = mean,geom="bar",color="black",width=0.55,size=0.3)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.08)+
  labs(x="Time(month)",y="Normalized stochastic ratio (%)")+
  scale_fill_manual(values=c("#d77c2b"))+
  theme_test()+theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
                     axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
                     axis.text.y=element_text(size=12,colour = 'black'),
                     panel.border = element_rect(size=1),
                     legend.text = element_text(size=10),
                     legend.title = element_text(size=10))+ylim(0,1)+
  theme(legend.position = "none")


#Figure 5b mantel tests
env<-read.csv('Environment.csv',row.names = 1) # sample as row, env as column
OTUs<-t(read.delim('eukspecies-time-normalised.txt',row.names = 1))
env<-apply(env[,1:10],2,as.numeric) #transfer numeric
OTUs<-vegdist(OTUs,method = "bray")
env<-vegdist(env,method = 'euclidean',na.rm = T) 
#partial mantel test
mantel(OTUs,env,permutations = 9999,method="spearman",na.rm=T)

df_mantel <- mantel_test(OTUs, env, spec.select = list("spc01" = 1:36))
df_mantel<- df_mantel %>% mutate(df_r=cut(r,breaks=c(-Inf,0.1,0.2,Inf),
                                          labels=c("<0.1","0.1-0.2",">=0.2")),
                                 df_p=cut(p.value,breaks=c(-Inf,0.001,0.05,Inf),
                                          labels=c("<0.001","0.001-0.05",">=0.05")))
df_mantel$linetype <- ifelse(df_mantel$p.value>=0.05,2,1)
df_mantel$linetype <- factor(df_mantel$linetype,levels = c("1","2"))
write.table(df_mantel,"df_mantel.txt", sep = '\t', col.names = NA, quote = FALSE)

quickcor(env, type = "lower",method = "spearman",show.diag = F,cor.test = T) +
  geom_square()+  
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"), limits = c(-1, 1),breaks = seq(-1,1,0.5))+
  geom_square() +
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r, 
                           linetype = linetype),
            label.size =4,
            label.fontface =1,
            curvature =0.2,
            nudge_x =0.2)+
  scale_size_manual(values = c(0.2, 0.4, 1))+
  scale_colour_manual(values=c("#87CEEB","gray","#9ACD32")) +
  scale_linetype_manual(values = c(1,2))+
  guides(fill = guide_colorbar(title = "Spearman's rho", order = 3),
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 1),
         linetype = "none") 