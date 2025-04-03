library(ggplot2)
library(dplyr)
library(cols4all)
library(patchwork)
library(grid)
library(eulerr)
library(vegan)
library(ggpubr)
library(ggsignif)


#Figure 1a
setwd('E:/time-meta')
df <- read.csv('euk-reads.csv',header = T)
head(df)

p <- ggplot(data = df, aes(x = time, y = counts, fill = Group)) +
  geom_bar(position = "dodge", stat="identity", color = NA, width = 0.8)
p

dt <- df %>%
  group_by(time, Group) %>% 
  summarise(m = mean(counts), 
            s = sd(counts)/sqrt(length(counts)))
head(dt)
dt$time <- factor(dt$time,levels = unique(dt$time))
dt$Group <- factor(dt$Group, levels = c('EukRep', 'Tiara', 'CCMetagen'))

p1 <- ggplot(dt, aes(x = time, y = m, fill = Group)) +
  geom_bar(position = position_dodge(), stat="identity", color = NA, width = 0.8) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  geom_errorbar(aes(ymin = m-s,
                    ymax = m+s),
                position = position_dodge(width = 0.8), 
                width = 0.2) +
  labs(y = 'Eukaryotic reads number')
p1


p2 <- p1 +
  scale_fill_manual(values = c("#b03525", "#1b6baa", "#d77c2b")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p2

p3 <- p2 + coord_cartesian(ylim = c(0,20000)) +
  theme_classic()+
  theme(legend.position = "none")
p3


p4 <- p2 + coord_cartesian(ylim = c(500000,2000000)) +
  theme_classic() +
  theme(axis.line.x = element_line(colour = "white"),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
  )
p4

grid.newpage()
plot_site1 <- viewport(x = 0.008, y = 0, width =0.9, height = 0.4, just = c('left','bottom'))
plot_site2 <- viewport(x = 0.008, y = 0.4, width =1.1, height = 0.4, just = c('left','bottom'))
print(p3,vp = plot_site1)
print(p4,vp = plot_site2)


#Figure 1b
venn_dat  <- read.delim("Venn.txt")
venn_list <- as.list(venn_dat)                         
venn_list <- purrr::map(venn_list, na.omit)            
venn_list <- lapply(venn_list, function(x) x[x != ""]) 
venn_list <- lapply(venn_list, unique)                 
venn_list

plot(euler(
  venn_list,
  shape = "circle"),                    
  quantities = list(type = c("percent","counts"),cex=1),          
  labels=list(cex=1),                   
  edges = list(col = "black", lex = 2), 
  fills = list(fill = c("#b03525","#1b6baa","#d77c2b"),alpha=0.7) 
)


#Figure 1c
data <-read.csv('euk-Kraken2.csv')
head(data)
dt <- data %>%
  group_by(time, Group) %>% 
  summarise(m = mean(counts), 
            s = sd(counts)/sqrt(length(counts)))
head(dt)
dt$Group <- factor(dt$Group, levels = c('EukRep', 'Tiara', 'CCMetagen'))
p <- ggplot(dt, aes(x = time, y = m, color = Group)) +
  geom_point()+
  geom_line(position = position_dodge(0.1),cex=1)+
  geom_errorbar(aes(ymin = m - s, ymax = m + s), 
                width = 0.2,cex=1)+
  scale_x_continuous(expand = c(0,0),limits = c(7,19))+
  scale_y_continuous(expand = c(0,0),limits = c(0,100))+
  scale_color_manual(values = c("#b03525","#1b6baa","#d77c2b"))+
  labs(x='Time (month)',y='Consistency ratio(%)')+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.8,0.55),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(size = 1,fill = 'transparent'))
p


#Figure 1d NMDS
library(ggplot2)
otu <- read.csv("*species-normalised.csv",header = TRUE ,row.names=1,sep = ",",stringsAsFactors = TRUE)#Species-normalised.csv of different approaches
otu <- t(otu)

otu.distance <- vegdist(otu, method = 'bray')
df_nmds <- metaMDS(otu.distance, k = 2)
summary(df_nmds)

otu <- read.csv("species-nmds.csv",header = TRUE ,row.names=1,sep = ",",stringsAsFactors = TRUE)#Combined table after calculating Bray-Curtis distances for the three methods
group <- read.table("group-3.txt", sep='\t', header=T)
colnames(group) <- c("samples","group")

##PERMANOVA
ad.result <- adonis2(otu ~ group, data = group, permutations = 999, distance = 'bray',by="margin")   
ad.result
df <- merge(df_points,group,by="samples")
head(df)

df_nmds_stress <- df_nmds$stress
df_nmds_stress
stressplot(df_nmds)

df_points <- as.data.frame(df_nmds$points)
df_points$samples <- row.names(df_points)
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)


p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3)+
  theme_bw()
p

df$group = factor(df$group, levels = c("EukRep","Tiara","CCMetagen"))
color=c("#b03525","#1b6baa","#d77c2b")
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+
  theme_bw()+
  geom_point(aes(color = group), shape = 19, size=2)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_text(aes(label=samples, y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F)+
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())+
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))
p1

p2 <- ggplot(df,aes(x=group,y=NMDS2))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#b03525","#1b6baa","#d77c2b"))+
  geom_signif(comparisons = list(c("EukRep","Tiara"),
                                 c("EukRep","CCMetagen"),
                                 c("Tiara","CCMetagen")),
              map_signif_level = T, 
              test = t.test,
              y_position = c(0.3,0.4,0.35),
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),
              size=0.8,color="black")
p2

p3 <- ggplot(df,aes(x=group,y=NMDS1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  coord_flip()+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#b03525","#1b6baa","#d77c2b"))+
  geom_signif(comparisons = list(c("EukRep","Tiara"),
                                 c("EukRep","CCMetagen"),
                                 c("Tiara","CCMetagen")),
              map_signif_level = T,
              test = t.test,
              y_position = c(0.48,0.62,0.55),
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),
              size=0.8,color="black")
p3

p4 <- ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,4), align = "hv")
p4