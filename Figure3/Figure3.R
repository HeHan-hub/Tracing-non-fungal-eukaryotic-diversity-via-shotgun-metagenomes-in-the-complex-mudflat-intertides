library(pcutils)
library(ggplot2)
library(cols4all)
library(vegan)
library(picante)
library(reshape)
library(ggpmisc)
library(phyloseq)
library(GUniFrac)
library(reshape2)
library(dplyr)
library(stats)
library(nlme)
library(ggpubr)
library(hillR)
library(readr)

#Figure 3a-b
otu <- read.csv("euksphylum-time-normalised.csv", header=T, row.names = 1)
#otu <- read.csv("eukgenus-time-normalised.csv", header=T, row.names = 1)

metadata = read.table("group.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
head(metadata, n = 3)
metadata$Group = factor(metadata$Group, levels = c("Aug","Oct","Dec","Feb","April","June"))

p <- stackplot(otu, metadata, group = "Group",relative = T, flow = T, others=F,topN= 12)+
  scale_fill_manual(name="otu", values=c('#8bc8c0','#ec7a70','#b0ce6c',"#c9e1c1",'#FFB6C1','#7B68EE','#B0C4DE','#778899','#7daac9','#B0E0E6','#bebcbf'))
p


#Figure 3d TTR
otu <- read.csv('eukspecies-time-normalised.csv', row.names = 1)
otu <- t(otu)
observed_species <- rowSums(otu > 0)
observed_species
write.csv(observed_species, 'Richness.csv', quote = FALSE)

#Richness totalization
richness <- read.csv('Richness.csv', row.names = 1)
time<- read.csv("time.csv", row.names = 1)
time <- as.data.frame(time)

aa <- cbind(time, richness)
aa <- log10(aa)
write.csv(aa, 'aa.csv', quote = FALSE)
aa <- read.csv('aa.csv',row.names = 1)
summary(lm(aa$XRichness ~ aa$time))
summary(lm(aa$Richness18 ~ aa$time))

data2=melt(data.frame(aa), id="time")
data2

model <- lm(value ~ time, data = data2)
confint(model, level = 0.95)  

p <- ggplot(data=data2, aes(x=time, y=value, group=variable, color=variable, shape= variable, size = variable)) +
  scale_color_manual("", values = c("#b03525", "#1b6baa"))+
  labs(x =expression("log(Time(month))"),y="log(Cumulative Richness)",
       title = "Taxa¨Ctime relationships")+
  scale_shape_manual(values = c(19,17))+
  scale_size_manual(values = c(2,2))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.justification=c(1,0), legend.position=c(1,0))+
  ylim(2,4)+
  geom_point()+
  geom_smooth(method = "lm", level=0.95)+
  stat_poly_eq(aes(x = time, y = value,label = paste(..rr.label..,..p.value.label.., sep = '~~~~')), formula = y ~ x, parse = T)+
  geom_text(aes(x=1,y=4,label="¦Ø="),
            color="black",family = "serif",fontface = "plain",size = 5)
p1


#Figure 3e TDR
time<- read.csv("time.csv", row.names = 1)
otu <- read.csv("eukspecies-time-normalised.csv", header=T, row.names = 1)

Bray.curtis <-as.matrix(vegdist(t(otu), method = 'bray'))
head(Bray.curtis)
similarity <- 1-Bray.curtis

similarity.reshape <- reshape2::melt(similarity)
colnames(similarity.reshape) <- c("Var1","Var2","similarity")

site_dis_df <- as.data.frame(time)
time_matrix <- as.matrix(site_dis_df)
distances <- vegdist(site_dis_df,method = "euclidean")
dist_matrix <- as.matrix(distances)
rownames(dist_matrix) <- rownames(time)
colnames(dist_matrix) <- rownames(time)
site_dis2 <- reshape2::melt(dist_matrix)
head(site_dis2)
colnames(site_dis2) <- c("Var1","Var2","distance")

final.table <- cbind(site_dis2,similarity.reshape)
plot.table <- as.data.frame(cbind(site_dis2$distance,similarity.reshape$similarity))
colnames(plot.table) <- c("distance","similarity")
plot.table <- log(plot.table)
write.csv(plot.table,"TDR.csv")
plot.table <- read.csv('TDR.csv', row.names = 1)

summary<-gls(similarity ~ poly(distance,2),data=plot.table,method="ML")
summary1<-gls(similarity ~ distance,data=plot.table)
summary
summary1
summary18<-gls(similarity18 ~ poly(distance,2),data=plot.table,method="ML")
summary181<-gls(similarity18 ~ distance,data=plot.table)
summary18
summary181

data2=melt(data.frame(plot.table), id="distance")
data2

show_point_shapes()
p <- ggplot(data=data2, aes(x=distance, y=value, group=variable, color=variable, shape= variable, size = variable)) +
  scale_color_manual("", values = c("#b03525", "#1b6baa"))+
  labs(x =expression("ln(Time(month))"),y="ln(1-Bray-curtis)",
       title = "Time¨Cdecay relationships")+
  scale_shape_manual(values = c(19,17))+
  scale_size_manual(values = c(3,2))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.justification=c(1,0), legend.position=c(1,0))+
  ylim(-1.2,0.2)+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_smooth(method = "lm", formula=y~poly(x,2))+
  stat_poly_eq(aes(x = distance, y = value,label = paste(..rr.label..,..p.value.label.., sep = '~~~~')), formula = y ~ x, parse = T)+
  geom_text(aes(x=2,y=0,label="slope="),
            color="black",family = "serif",fontface = "plain",size = 5)
p


#Figure 3f
comm <- read.csv("eukspecies-time-normalised.csv", row.names = 1)
comm <- t(comm)
q_values <- seq(0, 2, by = 0.1)
hill <- list()
for(q in q_values){
  hill[[as.character(q)]] <- hill_taxa_parti_pairwise(comm, q, rel_then_pool = FALSE, show_warning = TRUE, .progress = TRUE, output = "data.frame", pairs = "full")
}

write.csv(hill, "hill-b.csv")

data <- read.csv("beta-similarity-q.csv")#0¡Üq¡Ü2

values <- data[[2]]

n_rows <- 36
n_cols <- ceiling(length(values) / n_rows)  
filled_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)

for (i in 1:length(values)) {
  column_index <- (i - 1) %/% n_rows + 1  
  row_index <- (i - 1) %% n_rows + 1 
  filled_matrix[row_index, column_index] <- values[i]
}

filled_df <- as.data.frame(filled_matrix)

row_names <- data[[1]] 
col_names <- colnames(data)[-1] 

rownames(filled_df) <- row_names
colnames(filled_df) <- col_names[1:ncol(filled_df)]

print(filled_df)

write.csv(filled_df, "beta-similarity-q-1.csv", row.names = TRUE)

hill1 <- list() 
for(q in q_values){
  hill1[[as.character(q)]] <- hill_taxa(comm, q)
}

write.csv(hill1, "hill-a.csv")

data <-read.csv('hill-q.csv')
head(data)

p <- ggplot(data, aes(x = q)) +
  geom_point(aes(y = TTR), color = "#b03525") +
  geom_smooth(aes(y = TTR), method = "loess", color = "#b03525", se = FALSE) +  
  geom_point(aes(y = TDR), color = "#1b6baa") +  
  geom_smooth(aes(y = TDR), method = "loess", color = "#1b6baa", se = FALSE) +  
  labs(x = 'Diversity order(q)',y='Slope coefficient') +
  theme_minimal()
p