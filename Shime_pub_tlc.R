

#################################  PHYLOSEQ  #################################################################

##############################################################################################################
##############################################################################################################

setwd("C:/Users/latoe/Dropbox/PhD/Labmet/SHIME LNS/1.SHIME results/tryphyloseq/try2")

library(devtools)  # Load the devtools package

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(phyloseq)
library(knitr)
library(gtools)
library(scales)
library(ggthemes)
library("cowplot")
library(extrafont)
loadfonts(device = "win")

source("C:/Users/latoe/Dropbox/PhD/Labmet/SHIME LNS/1.SHIME results/tryphyloseq/try2/results_whole_treatment_period/Highqualgraphr.r")

# SET PLOTTING THEME
theme_set(theme_bw())
# Assign variables for imported data
sharedfile = "shime.shared"
taxfile = "shime.taxonomy"
mapfile = "MetaData.csv"
# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
# Import sample metadata
map0 <- read.csv(mapfile)
#convert the metadata to phyloseq sample data
map0$Treatment<-factor(map0$Treatment, labels=c("Control", "LNS", "LNSp"))
map0$Treatment.day <- with(map0, interaction(Treatment,  Day), drop = TRUE )
map0$idinfant<-factor(map0$idinfant, labels=c("Infant 1", "Infant 2", "Infant 3"))
map <- sample_data(map0)
# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID
# create a phyloseq object
moth_merge <- merge_phyloseq(mothur_data, map)
moth_merge
#assign taxonomic ranks
colnames(tax_table(moth_merge))
colnames(tax_table(moth_merge)) <- c("Domain", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

# Remove data from no treatment periods
Shime<- moth_merge %>%
  subset_samples(Period != "end_stabilization")
Shime<- Shime %>%
  subset_samples(Period != "wash_out" )


# Check phylogenic ranks
get_taxa_unique(Shime, "Phylum")
get_taxa_unique(Shime, "Class")
get_taxa_unique(Shime, "Order")
get_taxa_unique(Shime, "Family")
get_taxa_unique(Shime, "Genus")

### EXPLORE DATA ##

### CHECK DATA DISTRIBUTION ###

#' Make data frames with a column for the sum,
#' log sum, square root of sum and square sum 
#' of relative abundance of each sample
sample_sum_shime <- data.frame(sum = sample_sums(Shime))
sample_sum_shime$log10sum<-log10(sample_sum_shime$sum)
sample_sum_shime$sqrtsum<-sqrt(sample_sum_shime$sum)
sample_sum_shime$sqsum<-(sample_sum_shime$sum)^2

##' Density plots to check the shape of 
##' distribution of native and transformed data

hist1<-ggplot(sample_sum_shime, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "lightblue", binwidth = 500) +
  xlab("Relative abundance") +ylab("Frequency")+ 
  theme(panel.grid.major = element_blank(),text = element_text(size = 12), 
        panel.grid.minor = element_blank(),legend.title=element_blank())

# Histogram of sample log sums
hist2<-ggplot(sample_sum_shime, aes(x = log10sum)) + 
  geom_histogram(color = "black", fill = "lightblue", binwidth = 0.05) +
  xlab("Log sum relative abundance") + ylab("Frequency")+
  theme(panel.grid.major = element_blank(),text = element_text(size = 12), 
        panel.grid.minor = element_blank(),legend.title=element_blank())

# Histogram square root of sums
hist3<-ggplot(sample_sum_shime, aes(x = sqrtsum)) + 
  geom_histogram(color = "black", fill = "lightblue", binwidth = 5) +
  xlab("Square root sum of relative abundance ") + ylab("Frequency")+
  theme(panel.grid.major = element_blank(),text = element_text(size = 12), 
        panel.grid.minor = element_blank(),legend.title=element_blank())
hist3
# Histogram of square sums
hist4<-ggplot(sample_sum_shime, aes(x = sqsum)) + 
  geom_histogram(color = "black", fill = "lightblue", binwidth = 10000000) +
  xlab("square sum relative abundance") + ylab("Frequency")+
  theme(panel.grid.major = element_blank(),text = element_text(size = 12), 
        panel.grid.minor = element_blank(),legend.title=element_blank())
hist4
plot_grid(hist1, hist2, hist3, hist4, nrow = 2, ncol=2)
#### Data are skewed, logData are symetrical

# summary stats of sample relative abundances
smin <- min(sample_sums(Shime))
smean <- mean(sample_sums(Shime))
smax <- max(sample_sums(Shime))
smd <- median(sample_sums(Shime))
ssd<- sd(sample_sums(Shime))
## as expected, mean and median are different (skewed data)

#######################################################'
#####   Bar graphs of taxa      #######################
#######################################################
# melt to long format (for ggploting) 
# prune out phyla not detected in each sample
library("RColorBrewer")
Shime_phylum <- Shime %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.00) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


# Set colors for plotting
phylum_colors <- c("red","blue","orange","green","grey","purple", "yellow", 
                   "lightblue", "thistle2", "navy","pink", "red4", "black","lightgreen",  
                   "darkgreen",  "turquoise1", "violetred2"
)


# Plot 
barphylum<-ggplot(Shime_phylum, aes(x = Day, y = Abundance, fill = Phylum)) + 
  facet_grid(Treatment~idinfant) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 0%) \n") + 
  ggtitle("Phylum Composition Shime samples by treatment condition and per infant") 

barphylum+  theme(panel.grid = element_blank(),
                  text = element_text(family = "Arial", size = 12))

#plot Genera >1% abundance
Shime_Genus <- Shime %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                         
  arrange(Genus)

# plot of genus>1% wrapped per infant and treatment, with day and abundance as x and y axes)
bargenus<-ggplot(Shime_Genus, aes(x = Day, y = Abundance, fill = Genus)) + 
  facet_grid(Treatment~idinfant) +
  scale_x_continuous(breaks = round(seq(min(Shime_Genus$Day), max(Shime_Genus$Day), by = 3),1))+
  geom_bar(stat = "identity") + 
  xlab(label = "Day")+
  theme(panel.grid = element_blank(),text = element_text(family = "Arial", size = 16))+
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.y = element_text(vjust = 2), axis.text = element_text(colour = "black", size = 16))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus> 1%)") 
bargenus
highqualgraphR(bargenus,"figure_2",extension = "tiff",
               family="Arial")

### ANALYSIS OF ALPHA DIVERSITY ### this can be part of supplemental material
theme_set(theme_bw())

rich<-estimate_richness(Shime, split= TRUE, measures = NULL)

map1<-filter(map0, Period !="wash_out"& Period !="end_stabilization")

#map1$Treatment.day <- with(map1, interaction(treatment,  day), drop = TRUE )
rich$idinfant<-map1$idinfant
rich$day<-map1$Day
rich$treatment<-map1$Treatment
## Mean richness and evenness per treatment condition and time point(all infants together)
df1 <- rich %>% group_by(treatment, idinfant ) %>% summarise_at(vars(Observed, Shannon, InvSimpson, Chao1), sum)
df2 <- rich %>% group_by(treatment, day ) %>% summarise_at(vars(Observed, Shannon, InvSimpson, Chao1), sum)

Rich<-ggplot(df2, aes(x=day, y=Observed, col=treatment ))+geom_point(size=4)+ geom_line()+
  xlab("Day")+ ylab("Observed richness")+ 
  scale_x_continuous(breaks = round(seq(min(df2$day), max(df2$day), by = 3),1))+
  scale_color_manual(name="Treatment",values = c("red", "blue", "darkgreen"))+
  theme(panel.grid = element_blank(), legend.title = element_text("Treatment"),
        text = element_text(family = "Arial", size = 16),
        axis.title.y = element_text(vjust = 2), 
        axis.text = element_text(colour = "black", size = 16))

Shannon<-ggplot(df2, aes(x=day, y=Shannon, col=treatment ))+geom_point(size=4)+ geom_line()+
  xlab("Day")+ ylab("Shannon index")+ 
  scale_x_continuous(breaks = round(seq(min(df2$day), max(df2$day), by = 3),1))+
  scale_color_manual(name="Treatment", values = c("red", "blue", "darkgreen"))+
  theme(panel.grid = element_blank(), legend.title = element_text("Treatment"),
        text = element_text(family = "Arial", size = 16),
        axis.title.y = element_text(vjust = 2), 
        axis.text = element_text(colour = "black", size = 16))

invsimp<-ggplot(df2, aes(x=day, y=InvSimpson, col=treatment ))+geom_point(size=4)+ geom_line()+
  xlab("Day")+ ylab("Inverse Simpson index")+
  scale_x_continuous(breaks = round(seq(min(df2$day), max(df2$day), by = 3),1))+
  scale_color_manual(name="Treatment",values = c("red", "blue", "darkgreen"))+
  theme(panel.grid = element_blank(), legend.title =element_text("Treatment"),
        text = element_text(family = "Arial", size = 16),
        axis.title.y = element_text(vjust = 2), 
        axis.text = element_text(colour = "black", size = 16))

chao1<-ggplot(df2, aes(x=day, y=Chao1, col=treatment ))+geom_point(size=4)+ geom_line()+
  xlab("Day")+ ylab("Chao1 index")+ 
  scale_x_continuous(breaks = round(seq(min(df2$day), max(df2$day), by = 3),1))+
  scale_color_manual(name="Treatment",values = c("red", "blue", "darkgreen"))+
  theme(panel.grid = element_blank(), legend.title = element_text("Treatment"),
        text = element_text(family = "Arial", size = 16),
        axis.title.y = element_text(vjust = 2), 
        axis.text = element_text(colour = "black", size = 16))
highqualgraphR(Rich,"richness",extension = "tiff",
               family="Arial")
highqualgraphR(Shannon,"shannon",extension = "tiff",
               family="Arial")
highqualgraphR(invsimp,"invsimp",extension = "tiff",
               family="Arial")
highqualgraphR(chao1,"chao1",extension = "tiff",
               family="Arial")


Richplot<-grid.arrange(Rich, Shannon, invsimp, chao1, ncol=2, nrow = 2)
library(ggpubr)
plotrich<-as_ggplot(Richplot)
highqualgraphR(plotrich,"Diversity_plots",extension = "tiff",
               family="Arial")
# very instable, no clear trend in richness indexes

##Boxplots of diversity indexes

boxobs<-ggplot(rich, aes(x=treatment, y=Observed, fill=treatment ))+geom_boxplot()+
  ylab("Observed richness ")+scale_color_manual(values = c("red", "blue", "darkgreen")) +
  scale_fill_manual(values = c("red", "blue", "darkgreen"))+
  theme(text = element_text(family="Arial",size = 16, colour = "black"))+
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1))+
  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none")+
  scale_x_discrete("Treatment", labels=c("Control", "LNS", "LNSp"))+facet_wrap(~idinfant)
boxobs

boxsha<-ggplot(rich, aes(x=treatment, y=Shannon, fill=treatment ))+geom_boxplot()+
  ylab("Shannon Index")+scale_color_manual(values = c("red", "blue", "darkgreen")) +
  scale_fill_manual(values = c("red", "blue", "darkgreen"))+
  theme(text = element_text(family="Arial",size = 16, colour = "black"))+
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1))+
  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none")+
  scale_x_discrete("Treatment", labels=c("Control", "LNS", "LNSp"))+facet_wrap(~idinfant)

boxinvS<-ggplot(rich, aes(x=treatment, y=InvSimpson, fill=treatment ))+geom_boxplot()+
  ylab("Inverse Simpson Index")+scale_color_manual(values = c("red", "blue", "darkgreen")) +
  scale_fill_manual(values = c("red", "blue", "darkgreen"))+
  theme(text = element_text(family="Arial",size = 16, colour = "black"))+
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1))+
  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none")+
  scale_x_discrete("Treatment", labels=c("Control", "LNS", "LNSp"))+facet_wrap(~idinfant)

boxchao<-ggplot(rich, aes(x=treatment, y=Chao1, fill=treatment ))+geom_boxplot()+
  ylab("Chao1 Index")+scale_color_manual(values = c("red", "blue", "darkgreen")) +
  scale_fill_manual(values = c("red", "blue", "darkgreen"))+
  theme(text = element_text(family="Arial",size = 16, colour = "black"))+
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1))+
  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none")+
  scale_x_discrete("Treatment", labels=c("Control", "LNS", "LNSp"))+facet_wrap(~idinfant)
boxchao
highqualgraphR(boxobs,"boxrichness",extension = "tiff",
               family="Arial")
highqualgraphR(boxsha,"boxshannon",extension = "tiff",
               family="Arial")
highqualgraphR(boxinvS,"boxinvsimp",extension = "tiff",
               family="Arial")
highqualgraphR(boxchao,"boxchao1",extension = "tiff",
               family="Arial")
Richbox<-grid.arrange(boxobs, boxsha, boxinvS, boxchao, ncol=2, nrow = 2)
boxrich<-as_ggplot(Richbox)
highqualgraphR(boxrich,"Diversity_box",extension = "tiff",
               family="Arial")


###' boxes seem wider in the LNS treatment condition (maybe due to extremes values? ), 
###' however, confidence intervals overlap

######     STATISTICAL ANALYSIS OF COMMUNITY DIVERSITY     ######
#################################################################
library("lme4")
library("emmeans")
library("lmerTest")
library("lattice")

####BUILDING MODELS##########
### using blme to avoid singularity
library("blme")
# observed richness
mo0<-blmer(Observed~treatment+day +(1|idinfant),REML=FALSE, data=rich) 
summary(mo0)
mo01<-blmer(Observed~treatment+day +(1|treatment:idinfant),REML=FALSE, data=rich) 
summary(mo01)
mo02<-blmer(Observed~treatment+day + (treatment|idinfant),REML=FALSE, data=rich) #model is overfit

summary(mo02)
anova(mo0, mo01)
# model mo01 is best

#Shannon index
ms0<-blmer(Shannon~treatment+day +(1|idinfant),REML=FALSE, data=rich) 
summary(ms0)
ms01<-blmer(Shannon~treatment+day +(1|treatment:idinfant),REML=FALSE, data=rich)
summary(ms01)
ms02<-blmer(Shannon~treatment+day +(treatment|idinfant),REML=FALSE, data=rich) #model is overfit
summary(ms02)

anova(ms0, ms01,ms02)
#model ms01 is best
# inverse simpson index
mi0<-blmer(InvSimpson~treatment+day +(1|idinfant),REML=FALSE, data=rich) 
summary(mi0)
mi01<-blmer(InvSimpson~treatment+day +(1|treatment:idinfant),REML=FALSE, data=rich)
summary(mi01)
mi02<-blmer(InvSimpson~treatment+day +(treatment|idinfant),REML=FALSE, data=rich) #model is overfit
summary(mi02)
anova(mi0, mi01,mi02)
#model mi01 is best

#chao1 index
mc0<-blmer(Chao1~treatment+day +(1|idinfant),REML=FALSE, data=rich) 
summary(mc0)
mc01<-blmer(Chao1~treatment+day +(1|treatment:idinfant),REML=FALSE, data=rich)
summary(mc01)
mc02<-blmer(Chao1~treatment+day +(treatment|idinfant),REML=FALSE, data=rich)#model is overfit
summary(mc02)
anova(mc0, mc01,mc02)
# model mc0 and mc01 yield exactly the same results. 
msobs<-lmer(Observed~treatment+day +(1|treatment:idinfant),REML=TRUE, data=rich) 
mssha<-lmer(Shannon~treatment+day +(1|treatment:idinfant),REML=TRUE, data=rich) 
msinv<-lmer(InvSimpson~treatment+day +(1|treatment:idinfant),REML=TRUE, data=rich) 
mscha<-lmer(Chao1~treatment+day +(1|treatment:idinfant),REML=TRUE, data=rich)

summary(msobs)
summary(mssha)
summary(msinv)
summary(mscha)
### no difference in alpha diversity 


#######  ANALYSIS OF BETA DIVERSITY  #################
######################################################

######################################
#  Heat maps
####################################
Shime1<- Shime %>%subset_samples(idinfant=="Infant 1") 
Shime2<- Shime %>%subset_samples(idinfant=="Infant 2")
Shime3<- Shime %>%subset_samples(idinfant=="Infant 3")

##  20 most abundant taxa
Shime1_20 <- prune_taxa(names(sort(taxa_sums(Shime1),TRUE)[1:20]), Shime)

Shime2_20 <- prune_taxa(names(sort(taxa_sums(Shime2),TRUE)[1:20]), Shime2)

Shime3_20 <- prune_taxa(names(sort(taxa_sums(Shime3),TRUE)[1:20]), Shime3)

## Heatmap of the 20 most abundant taxa
dev.off()
(heatmap1<-plot_heatmap(Shime1_20, method = "NMDS", distance = "bray", 
                        sample.label = "Treatment.day",
                        taxa.label = "Genus" , low="blue", 
                        high="red", title = "Infant 1")+
    theme(panel.grid = element_blank(), 
          text = element_text(family = "Arial", size = 16),
          axis.title.y = element_text(vjust = 2), 
          axis.text = element_text(colour = "black", size = 16)))


(heatmap2<-plot_heatmap(Shime2_20, "NMDS", "bray", "Treatment.day", "Genus",   low="blue", 
                        high="red", title = "Infant 2")+
    theme(panel.grid = element_blank(), 
          text = element_text(family = "Arial", size = 16),
          axis.title.y = element_text(vjust = 2), 
          axis.text = element_text(colour = "black", size = 16)))


(heatmap3<-plot_heatmap(Shime3_20, "NMDS", "bray", "Treatment.day", "Genus",   low="blue", 
                        high="red", title = "Infant 3")+
    theme(panel.grid = element_blank(), 
          text = element_text(family = "Arial", size = 16),
          axis.title.y = element_text(vjust = 2), 
          axis.text = element_text(colour = "black", size = 16)))

heatmap<-ggarrange(heatmap1, heatmap2, heatmap3, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
heatmap
highqualgraphR(heatmap,"heatmaptab",extension = "tiff",
               family="Arial")
### the above graph tables are not readable for human eyes
### Making individual graphs per infant

highqualgraphR(heatmap1,"supplementary_figure 1",extension = "tiff",
               family="Arial")
highqualgraphR(heatmap2,"supplementary_figure 2",extension = "tiff",
               family="Arial")
highqualgraphR(heatmap3,"supplementary_figure 3",extension = "tiff",
               family="Arial")

######################################################
####  NMDS plots   ##############gt

# extract and format OTU tables 
abund_table1<-as.data.frame(t(Shime1@otu_table))
abund_table2<-as.data.frame(t(Shime2@otu_table))
abund_table3<-as.data.frame(t(Shime3@otu_table))

#format metadata per infant
meta_table1<-subset(map1, idinfant=="Infant 1")
meta_table2<-subset(map1, idinfant=="Infant 2")
meta_table3<-subset(map1, idinfant=="Infant 3")

row.names(meta_table1)=meta_table1$SampleID
row.names(meta_table2)=meta_table2$SampleID
row.names(meta_table3)=meta_table3$SampleID

meta_table1 = meta_table1[,-which(names(meta_table1) %in% c("SampleID"))]
meta_table2 = meta_table2[,-which(names(meta_table2) %in% c("SampleID"))]
meta_table3 = meta_table3[,-which(names(meta_table3) %in% c("SampleID"))]


###################################################################################################
shime.dist1<-metaMDS(abund_table1,distance = "bray", trymax = 100)
stressplot(shime.dist1)
shime.dist2<-metaMDS(abund_table2,distance = "bray", trymax = 100)
stressplot(shime.dist2)
shime.dist3<-metaMDS(abund_table3,distance = "bray", trymax = 100)
stressplot(shime.dist3)


#' Make a new data frame, with treatment, day and infant information there, 
#' to be used for coloring
grouping_info1<-data.frame(row.names = meta_table1$SampleName)
grouping_info1$Infant<-meta_table1$idinfant
grouping_info1$Treatment<-meta_table1$Treatment
grouping_info1$Day<-meta_table1$Day


grouping_info2<-data.frame(row.names = meta_table2$SampleName)
grouping_info2$Infant<-meta_table2$idinfant
grouping_info2$Treatment<-meta_table2$treatment
grouping_info2$Day<-meta_table2$day

grouping_info3<-data.frame(row.names = meta_table3$SampleName)
grouping_info3$Infant<-meta_table3$idinfant
grouping_info3$Treatment<-meta_table3$treatment
grouping_info3$Day<-meta_table3$day

head(grouping_info1)
head(grouping_info2)
head(grouping_info3)

Shime.nmds1=data.frame(x=shime.dist1$point[,1],y=shime.dist1$point[,2],Infant=as.factor(grouping_info1[,1]),Treatment=as.factor(grouping_info1[,2]),Sampling_Time=as.factor(grouping_info1[,3]))
Shime.nmds2=data.frame(x=shime.dist2$point[,1],y=shime.dist2$point[,2],Infant=as.factor(grouping_info2[,1]),Treatment=as.factor(grouping_info2[,2]),Sampling_Time=as.factor(grouping_info2[,3]))
Shime.nmds3=data.frame(x=shime.dist3$point[,1],y=shime.dist3$point[,2],Infant=as.factor(grouping_info3[,1]),Treatment=as.factor(grouping_info3[,2]),Sampling_Time=as.factor(grouping_info3[,3]))

graphics.off()
gg_envfit <- function(ord, env, groups=NA, scaling = 1, choices=c(1,2), perm = 999, alpha = 0.05, angle=20, len=0.5, unit="cm", arrow.col="red", pt.size=3, plot=TRUE) {
  df_ord <- vegan::scores(ord, display = "sites", choices = choices, scaling = scaling)
  df_ord <- as.data.frame(df_ord)
  axis.labels <- colnames(df_ord)
  if (!is.na(groups[1])) {
    df_ord$Group <- groups
    df_ord <- df_ord[ , c(3,1,2)]
    colnames(df_ord) <- c("Group", "x", "y")
  } else{
    colnames(df_ord) <- c("x", "y")
  }
  fit <- vegan::envfit(ord, env, choices = choices, perm = perm)
  if (min(fit$vectors$pvals) > alpha) {
    print(paste("No variable significant at alpha <=", as.character(alpha), sep = " "))
  } else {
    df_arrows <- as.data.frame(scores(fit, "vectors"))
    mult <- scale_arrow(df_arrows, df_ord[ , c("x", "y")])
    df_arrows <- mult * df_arrows
    df_arrows$var <- rownames(df_arrows)
    df_arrows$p.val <- fit$vectors$pvals
    colnames(df_arrows) <- c("x", "y", "var", "p.val")
    df_arrows <- df_arrows[df_arrows$p.val<=alpha, ]
    
    xlab <- axis.labels[1]
    ylab <- axis.labels[2]
    
    if (is.na(groups[1])) {
      plt <- ggplot(data=df_ord, aes(x=x, y=y)) + geom_point(size=pt.size) +
        xlab(xlab) + ylab(ylab)
    }
    else {
      plt <- ggplot(data=df_ord, aes(x=x, y=y, color=Group)) + geom_point(size=pt.size) +
        xlab(xlab) + ylab(ylab)
    }
    plt <- plt +
      geom_segment(data=df_arrows, aes(x=0, xend=x, y=0, yend=y),
                   arrow=arrow(angle=angle, length=unit(len, unit)), color=arrow.col) +
      geom_text(data=df_arrows, aes(x=x, y=y, label=var), color=arrow.col, hjust="outward")
    
    plt <- plt + coord_fixed(ratio=1)
    
    # Plot?
    if (plot) {print(plt)}
    
    # Return data frames, plot as a list.
    invisible(list(df_ord=df_ord, df_arrows=df_arrows, plot=plt))
  }
  
}
gg_ordibubble <- function(ord, env.var, var.label="Level", choices=c(1,2), plot=TRUE) {
  df_ord <- as.data.frame(vegan::scores(ord, display="sites", choices=choices))
  axis.labels <- colnames(df_ord)
  df_ord$var <- env.var
  xlab <- axis.labels[1]
  ylab <- axis.labels[2]
  colnames(df_ord) <- c("x", "y", var.label)
  
  plt <- ggplot(data=df_ord, aes(x=x, y=y, size=env.var)) +
    geom_point() + xlab(xlab) + ylab(ylab) + labs(size=var.label) +
    coord_fixed(ratio=1)
  
  # Plot?
  if (plot) {print(plt)}
  
  # Return data frames, plot as a list.
  invisible(list(df_ord=df_ord, plot=plt))
}

gg_ordicluster <- function (ord, cluster, treatments=NA, choices=c(1,2), prune = 0, col = 1, pt.size = 3, plot=TRUE)
{
  if (is.numeric(treatments)) {
    stop("'treatments' cannot be numeric")
  }
  n.trts <- nlevels(as.factor(treatments))
  if (n.trts==0) {
    treatments <- "none"
    n.trts <- 1
  }
  if (n.trts==1) {
    show.legend <- FALSE
  }
  else {
    show.legend <- TRUE
  }
  display <- "sites"
  w <- stats::weights(ord, display)
  weights.default <- function(object, ...) NULL
  w <- eval(w)
  mrg <- cluster$merge
  ord.scores <- scores(ord, display = display, choices=choices)
  if (nrow(mrg) != nrow(ord.scores) - 1)
    stop("Dimensions do not match in 'ord' and 'cluster'")
  if ((nrow(ord.scores) != length(treatments)) & (n.trts > 1))
    stop("Dimensions of 'ord' and 'treatments' do not match")
  if (length(w) == 1)
    w <- rep(w, nrow(ord.scores))
  n <- if (is.null(w))
    rep(1, nrow(ord.scores))
  else w
  noden <- numeric(nrow(mrg) - prune)
  go <- matrix(0, nrow(mrg) - prune, 2)
  col <- rep(col, length = nrow(ord.scores))
  col <- col2rgb(col)/255
  nodecol <- matrix(NA, nrow(mrg) - prune, 3)
  for (i in 1:(nrow(mrg) - prune)) {
    a <- mrg[i, 1]
    b <- mrg[i, 2]
    one <- if (a < 0)
      ord.scores[-a, ]
    else go[a, ]
    two <- if (b < 0)
      ord.scores[-b, ]
    else go[b, ]
    n1 <- if (a < 0)
      n[-a]
    else noden[a]
    n2 <- if (b < 0)
      n[-b]
    else noden[b]
    xm <- weighted.mean(c(one[1], two[1]), w = c(n1, n2))
    ym <- weighted.mean(c(one[2], two[2]), w = c(n1, n2))
    go[i, ] <- c(xm, ym)
    noden[i] <- n1 + n2
    colone <- if (a < 0)
      col[, -a]
    else nodecol[a, ]
    coltwo <- if (b < 0)
      col[, -b]
    else nodecol[b, ]
    nodecol[i, ] <- (n1 * colone + n2 * coltwo)/noden[i]
    
    # Rather than plotting the line segments, collect the coordinates and
    # color into a data frame.
    col.a = rgb(t(nodecol[i, ]))
    temp <- c(one[1], one[2], two[1], two[2], col.a)
    
    if (i==1){
      temp2 <- temp
    } else {
      temp2 <- rbind(temp2, temp)
      rownames(temp2) <- NULL # prevents duplicate row names
    }
    
  }
  
  colnames(temp2) <- c("x", "y", "xend", "yend", "Group")
  temp2 <- as.data.frame(temp2)
  j <- sapply(temp2, is.factor)
  temp2[j] <- lapply(temp2[j], as.character)
  j <- c( rep(TRUE, 4), FALSE)
  temp2[j] <- lapply(temp2[j], as.numeric)
  df_segments <- temp2
  
  df_ord <- as.data.frame(ord.scores)
  axis.labels <- colnames(df_ord)
  df_ord$Treatment <- treatments
  colnames(df_ord) <- c("x", "y", "Treatment")
  
  xlab <- axis.labels[1]
  ylab <- axis.labels[2]
  
  plt <- ggplot() +
    geom_segment(data=df_segments, aes(x=x, y=y, xend=xend, yend=yend),
                 color=df_segments$Group,
                 show.legend = FALSE) +
    geom_point(data=df_ord, aes(x=x, y=y, shape=Treatment), size=pt.size,
               show.legend = show.legend) +
    xlab(xlab) + ylab(ylab) + coord_fixed(ratio=1)
  
  if (plot==TRUE) {print(plt)}
  
  invisible(list(df_ord=df_ord, df_segments=df_segments, plot=plt))
}
gg_ordiplot <- function(ord, groups, scaling = 1, choices = c(1,2), kind = c("sd", "se", "ehull"), conf=NULL, show.groups="all", ellipse = TRUE, label = FALSE, hull = FALSE, spiders = FALSE, pt.size = 3, plot=TRUE) {
  groups <- as.factor(groups)
  if (show.groups[1]=="all") {
    show.groups <- as.vector(levels(groups))
  }
  
  # Get site coordinates to plot.
  df_ord <- vegan::scores(ord, display = "sites", scaling=scaling, choices=choices)
  axis.labels <- colnames(df_ord)
  df_ord <- data.frame(x=df_ord[ , 1], y=df_ord[ , 2], Group=groups)
  
  # Get ellipse centers to annotate.
  df_mean.ord <- aggregate(df_ord[,1:2], by=list(df_ord$Group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% show.groups, ]
  
  # Get parameters from the ordiellipse function.
  if (is.null(conf)) {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", label = label)
  } else {
    rslt <- vegan::ordiellipse(ord, groups=groups, display = "sites", scaling=scaling, choices=choices, kind = kind, show.groups = show.groups, draw = "none", conf = conf, label = label)
  }
  
  # Get points to plot for the ellipses.
  df_ellipse <- data.frame()
  for(g in show.groups) {
    df_ellipse <- rbind(df_ellipse, cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                                             vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, rslt[[g]]$scale))),Group=g))
  }
  colnames(df_ellipse) <- c("x", "y", "Group")
  df_ellipse <- df_ellipse[ , c(3,1,2)]
  
  # Make data frame for hulls.
  rslt.hull <- vegan::ordihull(ord, groups = groups, scaling = scaling, choices = choices, show.groups = show.groups, draw = "none")
  df_hull <- data.frame()
  df_temp <- data.frame()
  for (g in show.groups) {
    x <- rslt.hull[[g]][ , 1]
    y <- rslt.hull[[g]][ , 2]
    Group <- rep(g, length(x))
    df_temp <- data.frame(Group = Group, x=x, y=y)
    df_hull <- rbind(df_hull, df_temp)
  }
  
  # Make a data frame for the spiders.
  df_spiders <- df_ord
  df_spiders$cntr.x <- NA
  df_spiders$cntr.y <- NA
  for (g in show.groups) {
    df_spiders[which(df_spiders$Group==g), 4:5] <- df_mean.ord[which(df_mean.ord==g), 2:3]
  }
  df_spiders <- df_spiders[ , c(3,4,5,1,2)]
  df_spiders <- df_spiders[order(df_spiders$Group), ]
  df_spiders <- df_spiders[df_spiders$Group %in% show.groups, ]
  
  # Make basic ggplot with ellipses.
  xlab <- axis.labels[1]
  ylab <- axis.labels[2]
  plt <- ggplot2::ggplot() +
    geom_point(data=df_ord, aes(x=x, y=y, color=Group), size = pt.size) +
    xlab(xlab) + ylab(ylab)
  
  # Add ellipses.
  if (ellipse == TRUE) {
    plt <- plt + geom_path(data = df_ellipse, aes(x=x, y=y, color=Group), show.legend = FALSE)
  }
  
  # Add labels.
  if (label == TRUE) {
    plt <- plt + geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group, color=Group), show.legend = FALSE)
  }
  
  # Add hulls.
  if (hull == TRUE) {
    plt <- plt + geom_path(data=df_hull, aes(x=x, y=y, color=Group), show.legend = FALSE)
  }
  
  # Add spiders.
  if (spiders == TRUE) {
    plt <- plt + geom_segment(data=df_spiders, aes(x=cntr.x, xend=x, y=cntr.y, yend=y, color=Group), show.legend = FALSE)
  }
  
  plt <- plt + coord_fixed(ratio=1)
  
  # Plot?
  if (plot) {print(plt)}
  
  # Return data frames, plot as a list.
  invisible(list(df_ord=df_ord, df_mean.ord=df_mean.ord, df_ellipse=df_ellipse, df_hull=df_hull, df_spiders=df_spiders, plot=plt))
}
gg_ordisurf <- function(ord, env.var, choices=c(1,2), var.label="Level", binwidth, pt.size = 3, plot=TRUE) {
  # Extract ordisurf data for plotting
  ordi <- vegan::ordisurf(ord ~ env.var, plot=FALSE) #created the ordisurf object
  ordi.grid <- ordi$grid #extracts the ordisurf object
  ordi.data <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  ordi.data$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  df_surf <- data.frame(na.omit(ordi.data)) #gets rid of the nas
  
  # Extract site coordinates for plotting.
  df_ord <- as.data.frame(scores(ord, choices = choices, display = "sites"))
  axis.labels <- colnames(df_ord)
  colnames(df_ord) <- c("x", "y")
  
  # Make axis labels.
  xlab <- axis.labels[1]
  ylab <- axis.labels[2]
  
  # Calculate default binwidth
  if(missing(binwidth)) {
    r <- range(env.var)
    binwidth <- (r[2]-r[1])/15
  }
  
  ## Plotting in ggplot2
  plt <- ggplot(data=df_ord, aes(x=x, y=y)) + geom_point(size = pt.size) +
    
    xlab(xlab) + ylab(ylab) +
    stat_contour(data = df_surf, aes(x=x, y=y, z=z, color= ..level..), binwidth=binwidth) +
    labs(color=var.label) + coord_fixed(ratio=1)
  
  # Print plot
  if (plot) {print(plt)}
  
  # Return data frames, plot as a list.
  invisible(list(df_ord=df_ord, df_surf=df_surf, plot=plt))
}
scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  
  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}
library("ggplot2")
par('mar')
par(mar=c(1,1,1,1))
plot1<-gg_ordiplot(shime.dist1, meta_table1$Treatment,  scaling =1, choices = c(1,2), kind = "sd", 
                   conf=0.95, show.groups="all",
                   ellipse = TRUE,  label = FALSE, hull = FALSE, spiders = TRUE, pt.size = 3, plot=TRUE)

## extract the plot and edit it in ggplot2
plotinf1 <-plot1$plot
plotinf1 <-plotinf1 + theme_bw() +theme(panel.grid = element_blank(), 
                                        text = element_text(family = "Arial", size = 16),
                                        axis.title.y = element_text(vjust = 2), 
                                        axis.text = element_text(colour = "black", size = 16),
                                        panel.border=element_blank(),
                                        axis.line = element_line(color = "black"))+ 
  expand_limits(y = c(-0.8, 0.8), x = c(-0.8, 0.8))+
  scale_x_continuous(breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8)) +
  scale_y_continuous(breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8))+
  scale_color_manual(name="Treatment", values = c("red", "blue", "darkgreen"))+ggtitle("Infant 1")
plotinf1

plot2<-gg_ordiplot(shime.dist2, meta_table2$Treatment, scaling = 1,
                   choices = c(1,2), kind = "sd", 
                   conf=0.95, show.groups="all", 
                   ellipse = TRUE, label = FALSE, hull = FALSE, spiders = TRUE, pt.size = 3, plot=TRUE)
## extract the plot and edit it in ggplot2
plotinf2 <-plot2$plot
plotinf2 <-plotinf2 + theme_bw() +theme(panel.grid = element_blank(), 
                                        text = element_text(family = "Arial", size = 16),
                                        axis.title.y = element_text(vjust = 2), 
                                        axis.text = element_text(colour = "black", size = 16),
                                        panel.border=element_blank(),
                                        axis.line = element_line(color = "black"))+
  expand_limits(y = c(-0.8, 0.8), x = c(-0.8, 0.8))+
  scale_x_continuous(breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8)) +
  scale_y_continuous(breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8))+
  scale_color_manual(name="Treatment", values = c("red", "blue", "darkgreen"))+ggtitle("Infant 2")
plotinf2

plot3<-gg_ordiplot(shime.dist3, meta_table3$Treatment, scaling = 1, choices = c(1,2), kind = "sd", 
                   conf=0.95, show.groups="all", 
                   ellipse = TRUE, label = FALSE, hull = FALSE, spiders = TRUE, pt.size = 3, plot=TRUE)

## extract the plot and edit it in ggplot2
plotinf3 <-plot3$plot
plotinf3 <-plotinf3+ theme_bw() +theme(panel.grid = element_blank(), 
                                       text = element_text(family = "Arial", size = 16),
                                       axis.title.y = element_text(vjust = 2), 
                                       axis.text = element_text(colour = "black", size = 16),
                                       panel.border=element_blank(),
                                       axis.line = element_line(color = "black"))+  
  expand_limits(y = c(-0.8, 0.8), x = c(-0.8, 0.8))+
  scale_x_continuous(breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8)) +
  scale_y_continuous(breaks = c(-0.8, -0.4, 0.0, 0.4, 0.8))+
  scale_color_manual(name="Treatment", values = c("red", "blue", "darkgreen"))+ggtitle("Infant 3")
plotinf3
library("ggpubr")
plotnmds<-ggarrange(plotinf1, plotinf2, plotinf3, ncol=2, nrow = 2, common.legend = TRUE, legend = "right")

plotnmds
highqualgraphR(plotnmds,"nmds",extension = "tiff",
               family="Arial")


##############################################################################################################
##############################################################################################################

#################################     DESEQ2    ##############################################################

##############################################################################################################
##############################################################################################################

#' title: "Differential abundance testing for Laeticia"
#' author: "F.M. Kerckhof"



#+ setup, include = FALSE
#### load required packages ####
library(readxl)
library(dplyr)
library(DESeq2)
library(edgeR)
library(metagenomeSeq)
library(knitr)
library(ggplot2)
library(gtools)
library(scales)
library(extrafont)
library(ggthemes)
library("cowplot")
library("ggpubr")
font_import()
y
loadfonts(device = "win")
#### User defined functions ####
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#+ spacer

#' # Loading & preformatting the data

#### load the data ####
shared<-as.data.frame(Shime@otu_table)
meta<-filter(map0, Period !="wash_out"& Period !="end_stabilization")


row.names(meta_table1)=meta_table1$SampleID
row.names(meta_table2)=meta_table2$SampleID
row.names(meta_table3)=meta_table3$SampleID

meta_table1 = meta_table1[,-which(names(meta_table1) %in% c("SampleID"))]
meta_table2 = meta_table2[,-which(names(meta_table2) %in% c("SampleID"))]
meta_table3 = meta_table3[,-which(names(meta_table3) %in% c("SampleID"))]


#### preformat the data ###
sh.x <- as.data.frame(shared)
#remove absolute singletons
sh.x.ns <- sh.x[which(rowSums(sh.x)!=1),]

#### prevalence filtering ####
## prevalence = (fraction of samples in which an OTU is observed minimum 1 time)
minobs=1
prevalence <- apply(as.matrix(sh.x.ns),1,function(x,minobs){sum(x>=minobs)},minobs)/ncol(sh.x.ns)
prevalencefilter <- prevalence>0.05
sh.x.prevfilt <- sh.x.ns[prevalencefilter,]
##Read counts should exceed 0.5 times the number of samples
sh.x.prevfilt <- sh.x.prevfilt[rowSums(sh.x.prevfilt)>0.5*ncol(sh.x.prevfilt),]
sh.x.prevfilt.prop <- sweep(sh.x.prevfilt,2,colSums(sh.x.prevfilt),'/')
# deseq normalise ==> very dependent upon design
##DESeq by treatment

###generating factor 3 as interaction treatment/idinfant to be used for contrasts
meta$Treatment <- factor(meta$Treatment)
meta$idinfant <- factor(meta$idinfant)

meta$Factor3 <- with(meta, interaction(Treatment,  idinfant), drop = TRUE )

class(meta$Factor3)
class(meta$time)
meta$timef<-as.factor(meta$time)
class(meta$timef)#+ spacer
class(meta$Day)
meta$dayf<-factor(meta$Day)


#+ deseqformatting
deseqdata <- as.matrix(sh.x.prevfilt +1)


#####################################################################################"""""
##################################################################################""
#####################################################################################"""
#DESEQ WITH FACTOR3, correcting for treatment period
#######################################################################

deseqdata1 <- DESeqDataSetFromMatrix(deseqdata,colData=meta,design= ~ timef+Factor3)
deseqdata1 <- estimateSizeFactors(deseqdata1)
sizeFactors(deseqdata1)
sharedfiltereddeseq1 <- counts(deseqdata1,normalized=TRUE)
deseqdata1 <- estimateDispersions(deseqdata1,fitType="local")
testdiff1 <- DESeq(deseqdata1,test = "Wald",fitType = "local")


# comparing LNSp and control treatment groups then LNS and control treatment groups  per infant 
#using bonferoni correction
####################################################
#Infant 1
####################################################

#infant 1: LNSp vs control
testdiff.resLNSp1 <- results(testdiff1, contrast = c("Factor3", "LNSp.Infant 1", "Control.Infant 1"),  
                             alpha=0.016,pAdjustMethod = "BY")
testdiff.res.orderedLNSp1 <- testdiff.resLNSp1[order(testdiff.resLNSp1$padj),]
summary(testdiff.resLNSp1)
### Results: 7 OTUs with positive LFC and 10 with negative LFC, no outliers
signotuLNSp1 <- testdiff.resLNSp1$padj<0.016

resdfLNSp1 <- data.frame(sh.x.prevfilt,sigdif=signotuLNSp1,padjusted=testdiff.resLNSp1$padj)
write.csv2(resdfLNSp1,"DELNSpinf1.csv")

mcols(testdiff.res.orderedLNSp1, use.names = TRUE)
write.csv2(testdiff.res.orderedLNSp1, "RELNSp1.csv") 

#infant 1: LNS vs control
testdiff.resLNS1 <- results(testdiff1,alpha=0.016, contrast = c("Factor3", "LNS.Infant 1", "Control.Infant 1"),
                            pAdjustMethod = "BY")
testdiff.res.orderedLNS1 <- testdiff.resLNS1[order(testdiff.resLNS1$padj),]
summary(testdiff.resLNS1)
#10 OTUs with positive LFC and 8 with positive LFC 
signotuLNS1 <- testdiff.resLNS1$padj<0.016
resdfLNS1 <- data.frame(sh.x.prevfilt,sigdif=signotuLNS1,padjusted=testdiff.resLNS1$padj)
write.csv2(resdfLNS1,"DELNSinf1.csv")
mcols(testdiff.res.orderedLNS1, use.names = TRUE)
write.csv2(testdiff.res.orderedLNS1, "RELNS1.csv") 

#infant 1: LNSp vs LNS
testdiff.resLNSp01 <- results(testdiff1,alpha=0.016, contrast = c("Factor3", "LNSp.Infant 1", "LNS.Infant 1"),
                              pAdjustMethod = "BY")
testdiff.res.orderedLNSp01 <- testdiff.resLNSp01[order(testdiff.resLNSp01$padj),]
summary(testdiff.resLNSp01)
#4 OTUs with positive LFC and 8 with positive LFC 
signotuLNSp01 <- testdiff.resLNSp01$padj<0.016
resdfLNSp01 <- data.frame(sh.x.prevfilt,sigdif=signotuLNSp01,padjusted=testdiff.resLNSp01$padj)
write.csv2(resdfLNSp01,"DELNSpinf01.csv")
mcols(testdiff.res.orderedLNSp01, use.names = TRUE)
write.csv2(testdiff.res.orderedLNSp01, "RELNSp01.csv") 


###################################################c
#Infant2
####################################################
#infant 2: LNSp vs control
testdiff.resLNSp2 <- results(testdiff1, contrast = c("Factor3", "LNSp.Infant 2", "Control.Infant 2"),  
                             alpha=0.016,pAdjustMethod = "BY")
testdiff.res.orderedLNSp2 <- testdiff.resLNSp2[order(testdiff.resLNSp2$padj),]
summary(testdiff.resLNSp2)
### Results: 7 OTUs are upregulated and 8 are downregulated, no outliers
signotuLNSp2 <- testdiff.resLNSp2$padj<0.016
resdfLNSp2 <- data.frame(sh.x.prevfilt,sigdif=signotuLNSp2,padjusted=testdiff.resLNSp2$padj)
write.csv2(resdfLNSp2,"DELNSpinf2.csv")
mcols(testdiff.res.orderedLNSp2, use.names = TRUE)
write.csv2(testdiff.res.orderedLNSp2, "RELNSp2.csv") 
#infant 2: LNS vs Control
testdiff.resLNS2 <- results(testdiff1,alpha=0.016, contrast = c("Factor3", "LNS.Infant 2", "Control.Infant 2"),
                            pAdjustMethod = "BY")
testdiff.res.orderedLNS2 <- testdiff.resLNS2[order(testdiff.resLNS2$padj),]
summary(testdiff.resLNS2)
#4 OTUs are upregulated and 7 are downregulated, no outliers
signotuLNS2 <- testdiff.resLNS2$padj<0.016
resdfLNS2 <- data.frame(sh.x.prevfilt,sigdif=signotuLNS2,padjusted=testdiff.resLNS2$padj)
write.csv2(resdfLNS2,"DELNSinf2.csv")
mcols(testdiff.res.orderedLNS2, use.names = TRUE)
write.csv2(testdiff.res.orderedLNS2, "RELNS2.csv") 


#infant 2: LNSp vs LNS
testdiff.resLNSp02 <- results(testdiff1,alpha=0.016, contrast = c("Factor3", "LNSp.Infant 2", "LNS.Infant 2"),
                              pAdjustMethod = "BY")
testdiff.res.orderedLNSp02 <- testdiff.resLNSp02[order(testdiff.resLNSp02$padj),]
summary(testdiff.resLNSp02)
#3 OTUs are upregulated and 3 are downregulated, no outliers
signotuLNSp02 <- testdiff.resLNSp02$padj<0.016
resdfLNSp02 <- data.frame(sh.x.prevfilt,sigdif=signotuLNSp02,padjusted=testdiff.resLNSp02$padj)
write.csv2(resdfLNSp02,"DELNSpinf02.csv")
mcols(testdiff.res.orderedLNSp02, use.names = TRUE)
write.csv2(testdiff.res.orderedLNSp02, "RELNSp02.csv") 
####################################################
#Infant3
####################################################

#infant 3: LNSp vs. control
testdiff.resLNSp3 <- results(testdiff1, contrast = c("Factor3", "LNSp.Infant 3", "Control.Infant 3"),  
                             alpha=0.016,pAdjustMethod = "BY")
testdiff.res.orderedLNSp3 <- testdiff.resLNSp3[order(testdiff.resLNSp3$padj),]
summary(testdiff.resLNSp3)
### Results: 6 OTUs are upregulated and 8 are downregulated, no outliers
signotuLNSp3 <- testdiff.resLNSp3$padj<0.016
resdfLNSp3 <- data.frame(sh.x.prevfilt,sigdif=signotuLNSp3,padjusted=testdiff.resLNSp3$padj)
write.csv2(resdfLNSp3,"DELNSpinf3.csv")
mcols(testdiff.res.orderedLNSp3, use.names = TRUE)
write.csv2(testdiff.res.orderedLNSp3, "RELNSp3.csv") 

#infant 3: LNS vs. control
testdiff.resLNS3 <- results(testdiff1,alpha=0.016, contrast = c("Factor3", "LNS.Infant 3", "Control.Infant 3"),
                            pAdjustMethod = "BY")
testdiff.res.orderedLNS3 <-testdiff.resLNS3[order(testdiff.resLNS3$padj),]
summary(testdiff.resLNS3)
#5 OTUs are upregulated and 6 are downregulated, no outliers
signotuLNS3 <- testdiff.resLNS3$padj<0.016
resdfLNS3 <- data.frame(sh.x.prevfilt,sigdif=signotuLNS3,padjusted=testdiff.resLNS3$padj)
write.csv2(resdfLNS3,"DELNSinf3.csv")
mcols(testdiff.res.orderedLNS3, use.names = TRUE)
write.csv2(testdiff.res.orderedLNS3, "RELNS3.csv") 

#infant 3: LNSp vs LNS
testdiff.resLNSp03 <- results(testdiff1,alpha=0.016, contrast = c("Factor3", "LNSp.Infant 3", "LNS.Infant 3"),
                              pAdjustMethod = "BY")
testdiff.res.orderedLNSp03 <- testdiff.resLNSp03[order(testdiff.resLNSp03$padj),]
summary(testdiff.resLNSp03)
#3 OTUs are upregulated and 6 are downregulated, no outliers
signotuLNSp03 <- testdiff.resLNSp03$padj<0.016
resdfLNSp03 <- data.frame(sh.x.prevfilt,sigdif=signotuLNSp03,padjusted=testdiff.resLNSp03$padj)
write.csv2(resdfLNSp03,"DELNSpinf03.csv")
mcols(testdiff.res.orderedLNSp03, use.names = TRUE)
write.csv2(testdiff.res.orderedLNSp03, "RELNSp03.csv") 
##############################################################################################

### GRAPHS OF SIGNIFICANT DE OTUs WITH LOG2 FOLD CHANGE >=|2|########

##############################################################################################
# Graphs infant 1
##############################################################################################


#################################################################################################

##    Differentially expressed OTUs LNSp vs. control, infant 1

#################################################################################################

otuDELNSpinf1 <- c("Otu00012",	"Otu00003",	"Otu00002",	"Otu00010",
                   "Otu00008",	"Otu00011",	"Otu00038",	"Otu00016",	
                   "Otu00014",	"Otu00028",	"Otu00013")


saveotu.DELNSpinf1 <- data.frame()
for(i in otuDELNSpinf1){
  diffotu1 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu1 <- cbind(diffotu1,OTU=i)
  saveotu.DELNSpinf1 <- rbind(saveotu.DELNSpinf1,diffotu1)
}


saveotu.DELNSpinf1$Factor2 <- plyr::mapvalues(saveotu.DELNSpinf1$Factor3,
                                              from=levels(saveotu.DELNSpinf1$Factor3),
                                              to=1:9)

saveotu.DELNSpinf1$OTUm <- factor(saveotu.DELNSpinf1$OTU,levels=c("Otu00012",	"Otu00003",	"Otu00002",	"Otu00010",
                                                                  "Otu00008",	"Otu00011",	"Otu00038",	"Otu00016",	
                                                                  "Otu00014",	"Otu00028",	"Otu00013"))
#make dataframe with only the contrast of interest
LNSpinf1<-saveotu.DELNSpinf1 %>%
  filter(Factor3 %in% c("LNSp.Infant 1", "Control.Infant 1"))

LNSpinf1$treatment<-factor(LNSpinf1$Factor2, labels = c("Control", "LNSp"))

LNSpinf1$OTUr<-factor(LNSpinf1$OTUm, labels = c("OTU12",	"OTU3",	"OTU2",	"OTU10",
                                                "OTU8",	"OTU11",	"OTU38",	"OTU16",	
                                                "OTU14",	"OTU28",	"OTU13"))
LNSpinf1$OTUm<-factor(LNSpinf1$OTUm, labels = c("Bifidobacterium \n (OTU12)",	"Megasphaera \n (OTU3)", 
                                                "Mitsuokella \n (OTU2)",	"Dialister \n (OTU10)",	
                                                "Bacteroides \n (OTU8)", "Acidaminococcus \n (OTU11)",	
                                                "Mitsuokella \n (OTU38)",	"Enterococcus \n (OTU16)",	
                                                "Enterobacteriaceae \n (OTU14)",
                                                "Burkholderiales \n (OTU28)",	
                                                "Veillonella \n (OTU13)"))


#################################################################################################

##    Differentially expressed OTUs LNS vs. control, infant 1

#################################################################################################

otuDELNSinf1 <- c("Otu00006",	"Otu00003",	"Otu00002",	"Otu00015",
                  "Otu00008",	"Otu00010",	"Otu00011",	"Otu00032",	
                  "Otu00030",	"Otu00050","Otu00013")


saveotu.DELNSinf1 <- data.frame()
for(i in otuDELNSinf1){
  diffotu1 <- plotCounts(testdiff1,gene=i, intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu1 <- cbind(diffotu1,OTU=i)
  saveotu.DELNSinf1 <- rbind(saveotu.DELNSinf1,diffotu1)
}



saveotu.DELNSinf1$Factor2 <- plyr::mapvalues(saveotu.DELNSinf1$Factor3,
                                             from=levels(saveotu.DELNSinf1$Factor3),
                                             to=1:9)

saveotu.DELNSinf1$OTUm <- factor(saveotu.DELNSinf1$OTU,levels=c("Otu00006",	"Otu00003",	"Otu00002",	"Otu00015",
                                                                "Otu00008",	"Otu00010",	"Otu00011",	"Otu00032",	
                                                                "Otu00030",	"Otu00050","Otu00013"))
#make dataframe with only the contrast of interest
LNSinf1<-saveotu.DELNSinf1 %>%
  filter(Factor3 %in% c("LNS.Infant 1", "Control.Infant 1"))

LNSinf1$treatment<-factor(LNSinf1$Factor2, labels = c("Control", "LNS"))

LNSinf1$OTUr<-factor(LNSinf1$OTUm, labels = c("OTU6",	"OTU3",	"OTU2",	"OTU15",
                                              "OTU8",	"OTU10",	"OTU11",	"OTU32",	
                                              "OTU30",	"OTU50","	OTU13"))


LNSinf1$OTUm<-factor(LNSinf1$OTUm, labels = c("Prevotella \n (OTU6)","	Megasphaera \n (OTU3)",	
                                              "Mitsuokella \n (OTU2)","Acinetobacter \n (OTU15)",	
                                              "Bacteroides \n (OTU8)",	"Dialister \n (OTU10)",
                                              "Acidaminococcus \n (OTU11)","Veillonella \n (OTU32)",	
                                              "Allisonella \n (OTU30)","Lactobacillus \n (OTU50)",	
                                              "Veillonella \n (OTU13)"))

#########################################################

##    Differentially expressed OTUs LNSp vs. LNS, infant 1

#################################################################################################

otuDELNSpinf01 <- c("Otu00012",	"Otu00004",	"Otu00014",	"Otu00035",
                    "Otu00032",	"Otu00013",	"Otu00015","Otu00006")


saveotu.DELNSpinf01 <- data.frame()
for(i in otuDELNSpinf01){
  diffotu1 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu1 <- cbind(diffotu1,OTU=i)
  saveotu.DELNSpinf01 <- rbind(saveotu.DELNSpinf01,diffotu1)
}


saveotu.DELNSpinf01$Factor2 <- plyr::mapvalues(saveotu.DELNSpinf01$Factor3,
                                               from=levels(saveotu.DELNSpinf01$Factor3),
                                               to=1:9)


saveotu.DELNSpinf01$OTUm <- factor(saveotu.DELNSpinf01$OTU,levels=c("Otu00012",	"Otu00004",	"Otu00014",	"Otu00035",
                                                                    "Otu00032",	"Otu00013",	"Otu00015","Otu00006"))
#make dataframe with only the contrast of interest
LNSpinf01<-saveotu.DELNSpinf01 %>%
  filter(Factor3 %in% c("LNSp.Infant 1", "LNS.Infant 1"))

LNSpinf01$treatment<-factor(LNSpinf01$Factor2, labels = c("LNS", "LNSp"))

LNSpinf01$OTUr<-factor(LNSpinf01$OTUm, labels = c("OTU12","OTU4","OTU14","OTU35",
                                                  "OTU32","OTU13","OTU15","OTU6"))


LNSpinf01$OTUm<-factor(LNSpinf01$OTUm, labels = c("Bifidobacterium \n (OTU12)","Escherichia/ \n Shigella(OTU4)",
                                                  "Enterobacteriaceae -\n (OTU14)",	"Veillonella \n (OTU35)",	
                                                  "Veillonella \n (OTU32)",	"Veillonella \n (OTU13)","Acinetobacter \n (OTU15)",
                                                  "Prevotella \n (OTU6)"))

##############################################################################################
##############################################################################################

# Graphs infant 2

##############################################################################################


#################################################################################################

##    Differentially expressed OTUs LNSp vs. control, infant 2

#################################################################################################

otuDELNSpinf2 <- c("Otu00019",	"Otu00009",	"Otu00004",	"Otu00010",	
                   "Otu00034",	"Otu00022",	"Otu00007",	"Otu00011",	
                   "Otu00008",	"Otu00025",	"Otu00024",	"Otu00018")


saveotu.DELNSpinf2 <- data.frame()
for(i in otuDELNSpinf2){
  diffotu2 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu2 <- cbind(diffotu2,OTU=i)
  saveotu.DELNSpinf2 <- rbind(saveotu.DELNSpinf2,diffotu2)
}


saveotu.DELNSpinf2$Factor2 <- plyr::mapvalues(saveotu.DELNSpinf2$Factor3,
                                              from=levels(saveotu.DELNSpinf2$Factor3),
                                              to=1:9)


saveotu.DELNSpinf2$OTUm <- factor(saveotu.DELNSpinf2$OTU,levels=c("Otu00019",	"Otu00009",	"Otu00004",	"Otu00010",	
                                                                  "Otu00034",	"Otu00022",	"Otu00007",	"Otu00011",	
                                                                  "Otu00008",	"Otu00025",	"Otu00024",	"Otu00018"))
LNSpinf2<-saveotu.DELNSpinf2 %>%
  filter(Factor3 %in% c("LNSp.Infant 2", "Control.Infant 2"))

LNSpinf2$treatment<-factor(LNSpinf2$Factor2, labels = c("Control", "LNSp"))

LNSpinf2$OTUr<-factor(LNSpinf2$OTUm, labels = c("OTU19",	"OTU9",	"OTU4",	"OTU10",	
                                                "OTU34",	"OTU22",	"OTU7",	"OTU11",	
                                                "OTU8",	"OTU25",	"OTU24",	"OTU18"))

LNSpinf2$OTUm<-factor(LNSpinf2$OTUm, labels = c("Bifidobacterium \n (OTU19)","Pseudomonas \n (OTU9)","Escherichia/ \n Shigella(OTU4)",
                                                "Dialister \n (OTU10)",	"Desulfovibrio \n (OTU34)",	"Stenotrophomonas \n (OTU22)",
                                                "Megasphaera \n (OTU7)","Acidaminococcus \n (OTU11)",	"Bacteroides \n (OTU8)",	"Morganella \n (OTU25)",
                                                "Lactobacillus \n (OTU24)",	"Bacteroides \n (OTU18)"))




#################################################################################################

##    Differentially expressed OTUs LNS vs. control, infant 2

#################################################################################################

otuDELNSinf2 <- c("Otu00009",	"Otu00015",	"Otu00034",	"Otu00007",	
                  "Otu00008",	"Otu00024",	"Otu00025",	"Otu00022")


saveotu.DELNSinf2 <- data.frame()
for(i in otuDELNSinf2){
  diffotu2 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu2 <- cbind(diffotu2,OTU=i)
  saveotu.DELNSinf2 <- rbind(saveotu.DELNSinf2,diffotu2)
}

saveotu.DELNSinf2$Factor2 <- plyr::mapvalues(saveotu.DELNSinf2$Factor3,
                                             from=levels(saveotu.DELNSinf2$Factor3),
                                             to=1:9)


saveotu.DELNSinf2$OTUm <- factor(saveotu.DELNSinf2$OTU,levels=c("Otu00009",	"Otu00015",	"Otu00034",	"Otu00007",	
                                                                "Otu00008",	"Otu00024",	"Otu00025",	"Otu00022"))

LNSinf2<-saveotu.DELNSinf2 %>%
  filter(Factor3 %in% c("LNS.Infant 2", "Control.Infant 2"))

LNSinf2$treatment<-factor(LNSinf2$Factor2, labels = c("Control", "LNS"))

LNSinf2$OTUr<-factor(LNSinf2$OTUm, labels = c("OTU9",	"OTU15",	"OTU34",	"OTU7",	
                                              "OTU8",	"OTU24",	"OTU25",	"OTU22"))

LNSinf2$OTUm<-factor(LNSinf2$OTUm, labels = c("Pseudomonas \n (OTU9)",	"Acinetobacter \n (OTU15)",	
                                              "Desulfovibrio \n (OTU34)",	"Megasphaera \n (OTU7)",	
                                              "Bacteroides \n (OTU8)",	"Lactobacillus \n (OTU24)",	
                                              "Morganella \n (OTU25)",	"Stenotrophomonas \n (OTU22)"))


#################################################################################################

##    Differentially expressed OTUs LNSp vs. LNS, infant 2

#################################################################################################

otuDELNSpinf02 <- c("Otu00019","Otu00022","Otu00010",		
                    "Otu00018",	"Otu00015")


saveotu.DELNSpinf02 <- data.frame()
for(i in otuDELNSpinf02){
  diffotu1 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu1 <- cbind(diffotu1,OTU=i)
  saveotu.DELNSpinf02 <- rbind(saveotu.DELNSpinf02,diffotu1)
}


saveotu.DELNSpinf02$Factor2 <- plyr::mapvalues(saveotu.DELNSpinf02$Factor3,
                                               from=levels(saveotu.DELNSpinf02$Factor3),
                                               to=1:9)



saveotu.DELNSpinf02$OTUm <- factor(saveotu.DELNSpinf02$OTU,levels=c("Otu00019","Otu00022","Otu00010",		
                                                                    "Otu00018",	"Otu00015"))
#make dataframe with only the contrast of interest
LNSpinf02<-saveotu.DELNSpinf02 %>%
  filter(Factor3 %in% c("LNSp.Infant 2", "LNS.Infant 2"))

LNSpinf02$treatment<-factor(LNSpinf02$Factor2, labels = c("LNS", "LNSp"))

LNSpinf02$OTUr<-factor(LNSpinf02$OTUm, labels = c("OTU19","OTU22",	"OTU10",		
                                                  "OTU18","OTU15"	))

LNSpinf02$OTUm<-factor(LNSpinf02$OTUm, labels = c("Bifidobacterium \n (OTU19)","Stenotrophomonas \n (OTU22)"	,	
                                                  "Dialister \n (OTU10)","Bacteroides \n (OTU18)","Acinetobacter \n (OTU15)"))


##############################################################################################
#                        Graphs infant 3
##############################################################################################


#################################################################################################

##    Differentially expressed OTUs LNSp vs. control, infant 3

#################################################################################################

otuDELNSpinf3 <- c("Otu00021",	"Otu00012",	"Otu00009",	"Otu00008",	
                   "Otu00011",	"Otu00024",	"Otu00034",	"Otu00025",	
                   "Otu00020",	"Otu00013",	"Otu00018",	"Otu00017")

saveotu.DELNSpinf3 <- data.frame()
for(i in otuDELNSpinf3){
  diffotu3 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu3 <- cbind(diffotu3,OTU=i)
  saveotu.DELNSpinf3 <- rbind(saveotu.DELNSpinf3,diffotu3)
}


saveotu.DELNSpinf3$Factor2 <- plyr::mapvalues(saveotu.DELNSpinf3$Factor3,
                                              from=levels(saveotu.DELNSpinf3$Factor3),
                                              to=1:9)


saveotu.DELNSpinf3$OTUm <- factor(saveotu.DELNSpinf3$OTU,levels=c("Otu00021",	"Otu00012",	"Otu00009",	"Otu00008",	
                                                                  "Otu00011",	"Otu00024",	"Otu00034",	"Otu00025",	
                                                                  "Otu00020",	"Otu00013",	"Otu00018",	"Otu00017"))
LNSpinf3<-saveotu.DELNSpinf3 %>%
  filter(Factor3 %in% c("LNSp.Infant 3", "Control.Infant 3"))

LNSpinf3$treatment<-factor(LNSpinf3$Factor2, labels = c("Control", "LNSp"))

LNSpinf3$OTUr<-factor(LNSpinf3$OTUm, labels = c("OTU21","OTU12","OTU9",	"OTU8",	
                                                "OTU11","OTU24","OTU34","OTU25",
                                                "OTU20","OTU13","OTU18","OTU17"))


LNSpinf3$OTUm<-factor(LNSpinf3$OTUm, labels = c("Megasphaera \n (OTU21)","Bifidobacterium \n (OTU12)",
                                                "Pseudomonas \n (OTU9)",	"Bacteroides \n (OTU8)",	
                                                "Acidaminococcus \n (OTU11)","Lactobacillus \n (OTU24)",
                                                "Desulfovibrio \n (OTU34)","Morganella \n (OTU25)",
                                                "Bacteroides \n (OTU20)","Veillonella \n (OTU13)",
                                                "Bacteroides \n (OTU18)","Bilophila \n (OTU17)"))


#################################################################################################

##    Differentially expressed OTUs LNS vs. control, infant 3

#################################################################################################
otuDELNSinf3 <- c("Otu00021",	"Otu00015",	"Otu00011",	"Otu00009",	
                  "Otu00005",	"Otu00024",	"Otu00025",	"Otu00017")


saveotu.DELNSinf3 <- data.frame()
for(i in otuDELNSinf3){
  diffotu3 <- plotCounts(testdiff1,gene=i,intgroup =c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu3 <- cbind(diffotu3,OTU=i)
  saveotu.DELNSinf3 <- rbind(saveotu.DELNSinf3,diffotu3)
}

saveotu.DELNSinf3$Factor2 <- plyr::mapvalues(saveotu.DELNSinf3$Factor3,
                                             from=levels(saveotu.DELNSinf3$Factor3),
                                             to=1:9)


saveotu.DELNSinf3$OTUm <- factor(saveotu.DELNSinf3$OTU,levels=c("Otu00021",	"Otu00015",	"Otu00011",	"Otu00009",	
                                                                "Otu00005",	"Otu00024",	"Otu00025",	"Otu00017"))

LNSinf3<-saveotu.DELNSinf3 %>%
  filter(Factor3 %in% c("LNS.Infant 3", "Control.Infant 3"))

LNSinf3$treatment<-factor(LNSinf3$Factor2, labels = c("Control", "LNS"))

LNSinf3$OTUr<-factor(LNSinf3$OTUm, labels = c("OTU21","OTU15","OTU11","OTU9",
                                              "OTU5","OTU24","OTU25",	"OTU17"))

LNSinf3$OTUm<-factor(LNSinf3$OTUm, labels = c("Megasphaera \n (OTU21)","Acinetobacter \n (OTU15)",
                                              "Acidaminococcus \n (OTU11)","Pseudomonas \n (OTU9)",
                                              "Enterobacteriaceae \n (OTU5)","Lactobacillus \n (OTU24)",
                                              "Morganella \n (OTU25)",	"Bilophila \n (OTU17)"))



#################################################################################################

##    Differentially expressed OTUs LNSp vs. LNS, infant 3

#################################################################################################

otuDELNSpinf03 <- c( "Otu00012",	"Otu00008",	"Otu00018",	"Otu00014",	
                     "Otu00017",	"Otu00020",	"Otu00015")


saveotu.DELNSpinf03 <- data.frame()
for(i in otuDELNSpinf03){
  diffotu1 <- plotCounts(testdiff1,gene=i,intgroup = c("Day","Factor3"),returnData = TRUE,transform=FALSE)
  diffotu1 <- cbind(diffotu1,OTU=i)
  saveotu.DELNSpinf03 <- rbind(saveotu.DELNSpinf03,diffotu1)
}


saveotu.DELNSpinf03$Factor2 <- plyr::mapvalues(saveotu.DELNSpinf03$Factor3,
                                               from=levels(saveotu.DELNSpinf03$Factor3),
                                               to=1:9)


saveotu.DELNSpinf03$OTUm <- factor(saveotu.DELNSpinf03$OTU,levels=c("Otu00012",	"Otu00008",	"Otu00018",	"Otu00014",	
                                                                    "Otu00017",	"Otu00020",	"Otu00015"))
#make dataframe with only the contrast of interest
LNSpinf03<-saveotu.DELNSpinf03 %>%
  filter(Factor3 %in% c("LNSp.Infant 3", "LNS.Infant 3"))

LNSpinf03$treatment<-factor(LNSpinf03$Factor2, labels = c("LNS", "LNSp"))

LNSpinf03$OTUr<-factor(LNSpinf03$OTUm, labels = c("OTU12",	"OTU8",	"OTU18",	
                                                  "OTU14","OTU17", "OTU20",	
                                                  "OTU15"))

LNSpinf03$OTUm<-factor(LNSpinf03$OTUm, labels = c("Bifidobacterium \n (OTU12)",	"Bacteroides \n (OTU8)",
                                                  "Bacteroides \n (OTU18)",	"Enterobacteriaceae \n (OTU14)",
                                                  "Bilophila \n (OTU17)", "Bacteroides \n (OTU20)",	
                                                  "Acinetobacter \n (OTU15)"))

###################################

####  Plots DESEQ  #####
#################################


### plots infant1
plotlns1<-ggplot(LNSinf1, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSinf1$Day), max(LNSinf1$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1, Control vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))

plotlnsp1<-ggplot(LNSpinf1, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSpinf1$Day), max(LNSpinf1$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1, Control vs. LNSp")+ theme(plot.title = element_text(family = "Arial", size=16))


plotlnsp01<-ggplot(LNSpinf01, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSpinf01$Day), max(LNSpinf01$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1, LNSp vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))

### plots infant2
plotlns2<-ggplot(LNSinf2, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSinf2$Day), max(LNSinf2$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2, Control vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))

plotlnsp2<-ggplot(LNSpinf2, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSpinf2$Day), max(LNSpinf2$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2, LControl vs. LNSp")+ theme(plot.title = element_text(family = "Arial", size=16))


plotlnsp02<-ggplot(LNSpinf02, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSpinf02$Day), max(LNSpinf02$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2, LNSp vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))

### plots infant3
plotlns3<-ggplot(LNSinf3, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSinf3$Day), max(LNSinf3$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant3, Control vs. LNS")+
  theme(plot.title = element_text(family = "Arial", size=16))
plotlns3
plotlnsp3<-ggplot(LNSpinf3, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSpinf3$Day), max(LNSpinf3$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant3, Control vs. LNSp ")+
  theme(plot.title = element_text(family = "Arial", size=16))


plotlnsp03<-ggplot(LNSpinf03, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_x_continuous(breaks = round(seq(min(LNSpinf03$Day), max(LNSpinf03$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant3, LNSp vs. LNS")+
  theme(plot.title = element_text(family = "Arial", size=16))
### Plots are too large and hardly readable on a word or pdf document 


###  Plots for main article text
Otu1<-subset(LNSinf1, OTU =="Otu00006" | OTU=="Otu00013")
Otu2<-subset(LNSpinf1, OTU =="Otu00012" | OTU=="Otu00013")
Otu3<-subset(LNSpinf01, OTU =="Otu00012" | OTU=="Otu00013")

Otu4<-subset(LNSinf2, OTU =="Otu00009" |  OTU=="Otu00024")
Otu5<-subset(LNSpinf2, OTU =="Otu00019" |  OTU=="Otu00018")
Otu6<-subset(LNSpinf02, OTU =="Otu00019" | OTU=="Otu00018")

Otu7<-subset(LNSinf3, OTU =="Otu00021" |  OTU=="Otu00017")
Otu8<-subset(LNSpinf3, OTU =="Otu00021" | OTU=="Otu00017")
Otu9<-subset(LNSpinf03, OTU =="Otu00012" | OTU=="Otu00017")



plototu1<-ggplot(Otu1, aes(x=Day, y=count, colour=treatment)) +  geom_point(size=4) + 
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+ 
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu1$Day), max(Otu1$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1, LNS vs. Control")+ theme(plot.title = element_text(family = "Arial", size=16))

plototu2<-ggplot(Otu2, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) + 
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+
  theme(axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu2$Day), max(Otu2$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1, LNSp vs. Control")+ theme(plot.title = element_text(family = "Arial", size=16))


plototu3<-ggplot(Otu3, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu3$Day), max(Otu3$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1, LNSp vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))


plototu4<-ggplot(Otu4, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) + 
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial", size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+ 
  theme(axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu4$Day), max(Otu4$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2, LNS vs. Control")+ theme(plot.title = element_text(family = "Arial", size=16))


plototu5<-ggplot(Otu5, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) + 
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+ 
  theme(axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu5$Day), max(Otu5$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2, LNSp vs. Control")+ theme(plot.title = element_text(family = "Arial", size=16))


plototu6<-ggplot(Otu6, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) + 
  geom_line()+ scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+ 
  theme(axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu6$Day), max(Otu6$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2, LNSp vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))


### Infant3 ###
plototu7<-ggplot(Otu7, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) + 
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "blue"))+ 
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(),  legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu7$Day), max(Otu7$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 3, LNS vs. Contol")+ theme(plot.title = element_text(family = "Arial", size=16))

plototu8<-ggplot(Otu8, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) +
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family="Arial", size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "darkgreen"))+ 
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu8$Day), max(Otu8$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 3, LNSp vs. Control")+ theme(plot.title = element_text(family = "Arial", size=16))

plototu9<-ggplot(Otu9, aes(x=Day, y=count, colour=treatment)) + geom_point(size=4) + 
  geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("")+ xlab("Day")+ theme(text = element_text(family = "Arial", size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+ 
  theme( axis.text = element_text(colour = "black", size=16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_x_continuous(breaks = round(seq(min(Otu9$Day), max(Otu9$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 3, LNSp vs. LNS")+ theme(plot.title = element_text(family = "Arial", size=16))


Otuplot<- grid.arrange( plototu1, plototu4, plototu7,plototu2, plototu5, 
                        plototu8,plototu3, plototu6, plototu9,
                        ncol= 3, nrow=3)

p1<-as_ggplot(Otuplot)
p1
highqualgraphR(p1, "figure 4",extension = "tiff",
               family="Arial")

### supplementary DESEq plots

OtuS1<-subset(LNSinf1, OTU !="Otu00006" & OTU!="Otu00013")
OtuS2<-subset(LNSpinf1, OTU !="Otu00012" & OTU!="Otu00013")
OtuS3<-subset(LNSpinf01, OTU !="Otu00012" & OTU!="Otu00013")

OtuS4<-subset(LNSinf2, OTU !="Otu00009" &  OTU!="Otu00024")
OtuS5<-subset(LNSpinf2, OTU !="Otu00019" &  OTU!="Otu00018")
OtuS6<-subset(LNSpinf02, OTU !="Otu00019" & OTU!="Otu00018")

OtuS7<-subset(LNSinf3, OTU !="Otu00021" &  OTU!="Otu00017")
OtuS8<-subset(LNSpinf3, OTU !="Otu00021" & OTU!="Otu00017")
OtuS9<-subset(LNSpinf03, OTU !="Otu00012" & OTU!="Otu00017")

plotLNSpinf1.wrap<-ggplot(OtuS2, aes(x=Day, y=count, colour=treatment,label_value(OTUm, multi_line = TRUE), group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+ 
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS2$Day), max(OtuS2$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1")

plotLNSinf1.wrap<-ggplot(OtuS1, aes(x=Day, y=count, colour=treatment, group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+
  theme(axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS1$Day), max(OtuS1$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1")


plotLNSpinf01.wrap<-ggplot(OtuS3, aes(x=Day, y=count, colour=treatment,label_value(OTUm, multi_line = TRUE), group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS3$Day), max(OtuS3$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 1")


plotLNSpinf2.wrap<-ggplot(OtuS5, aes(x=Day, y=count, colour=treatment,label_value(OTUm, multi_line = TRUE), group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial", size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+ 
  theme(axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS5$Day), max(OtuS5$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2")


plotLNSinf2.wrap<-ggplot(OtuS4, aes(x=Day, y=count, colour=treatment, group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+ 
  theme(axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS4$Day), max(OtuS4$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2")


plotLNSpinf02.wrap<-ggplot(OtuS6, aes(x=Day, y=count, colour=treatment,label_value(OTUm, multi_line = TRUE), group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x, n=4))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+ 
  theme(axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS6$Day), max(OtuS6$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 2")


### Infant3 ###
plotLNSpinf3.wrap<-ggplot(OtuS8, aes(x=Day, y=count, colour=treatment,label_value(OTUm, multi_line = TRUE), group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial",size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNSp"), values=c("red", "darkgreen"))+ 
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS8$Day), max(OtuS8$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 3")


plotLNSinf3.wrap<-ggplot(OtuS7, aes(x=Day, y=count, colour=treatment, group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family="Arial", size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("Control", "LNS"), values=c("red", "blue"))+ 
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS7$Day), max(OtuS7$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 3")

plotLNSpinf03.wrap<-ggplot(OtuS9, aes(x=Day, y=count, colour=treatment,label_value(OTUm, multi_line = TRUE), group=interaction(OTUm, treatment))) + 
  geom_point(size=4) + geom_line()+scale_y_continuous(labels=scientific,trans='log10',breaks = trans_breaks('log10',function(x) 10^x))+ 
  ylab("Normalized log transformed counts")+ xlab("Day")+ 
  theme(text = element_text(family = "Arial", size=16))+ 
  scale_color_manual(name= "Treatment",breaks=c("LNS", "LNSp"), values=c("blue", "darkgreen"))+ 
  theme( axis.text = element_text(colour = "black", size = 16), axis.line = element_line(color="black", size = 1.5))+
  theme(panel.grid = element_blank(), axis.title.y = element_text(vjust = 4))+
  scale_x_continuous(breaks = round(seq(min(OtuS9$Day), max(OtuS9$Day), by = 3),1))+
  facet_wrap(~OTUm,scales="free_x")+ggtitle("Infant 3")
library("gridExtra")
library("gtable")

DEplotLNSp<-grid.arrange(plotLNSpinf1.wrap, plotLNSpinf2.wrap, plotLNSpinf3.wrap, nrow=2, ncol=2)

DEplotLNS<-grid.arrange(plotLNSinf1.wrap, plotLNSinf2.wrap, plotLNSinf3.wrap, nrow = 2, ncol = 2)

DEplotLNSp1<-grid.arrange(plotLNSpinf01.wrap, plotLNSpinf02.wrap, plotLNSpinf03.wrap, ncol = 2, nrow = 2)

p2<-as_ggplot(DEplotLNSp)
p3<-as_ggplot(DEplotLNS)
p4<-as_ggplot(DEplotLNSp1)
p2


highqualgraphR(p2,"supplementary figure 4",extension = "tiff",
               family="Arial")
highqualgraphR(p3,"supplementary figure 5",extension = "tiff",
               family="Arial")
highqualgraphR(p4,"supplementary figure 6",extension = "tiff",
               family="Arial")

citation("phyloseq")
citation("vegan")
citation("pairwiseAdonis")
citation("DESeq2")
####################################################################################################
##                                    qPCR data analysis
####################################################################################################

library("dplyr")
library("tidyr")
library("ggpubr")
library("car")
library("plyr")
library(ggplot2)
# load data and take the mean of replicates
pcrrep<-read.csv('qpcrdata.csv', header = TRUE, sep = ",")
pcr<- ddply(pcrrep, c("SampleID", "idinfant", "treatment", "day"), summarise,
            Nrep    = sum(!is.na(Quantity)),
            meanrep = mean(Quantity, na.rm=TRUE),
            sdrep   = sd(Quantity, na.rm=TRUE),
            serep   = sdrep / sqrt(Nrep) 
)

## explore data and format variables
theme_set(theme_bw())
class(pcr$day)
class(pcr$treatment)

class(pcr$idinfant)
class(pcr$meanrep)
class(pcr$sdrep)
class(pcr$serep)

pcr$idinf<-factor(pcr$idinfant, labels=c("Infant 1", "Infant 2", "Infant 3"))
pcr$treat<-factor(pcr$treatment, labels=c("Control", "LNS", "LNSp"))
pcr$Period<-factor(pcr$day, labels=c("Day0", "Day14", "Day15", "Day21", "Day27"))
class(pcr$idinf)
class(pcr$treat)
class(pcr$Period)
###SUMMARY STATISTICS###
#Summary statistics number of gene copies in total, then by treatment condition
dplyr::tbl_df(pcr)

pcrtreat<- ddply(pcr, c("treatment"), summarise,
                 N    = sum(!is.na(meanrep)),
                 mean = mean(meanrep, na.rm=TRUE),
                 median=median(meanrep),
                 max=max(meanrep),
                 min= min(meanrep),
                 sd   = sd(meanrep, na.rm=TRUE),
                 se   = sd / sqrt(N)
)
pcrtreat$mean<- formatC(pcrtreat$mean,  format = "e", digits = 2)
pcrtreat$sd<-formatC( pcrtreat$sd,  format = "e", digits = 2)
pcrtreat$se<-formatC( pcrtreat$se,  format = "e", digits = 2)
pcrtreat$median<- formatC(pcrtreat$median,  format = "e", digits = 2)
pcrtreat$max<-formatC( pcrtreat$max,  format = "e", digits = 2)
pcrtreat$min<-formatC( pcrtreat$min,  format = "e", digits = 2)

pcrinf<- ddply(pcr, c("treatment", "idinfant"), summarise,
               N    = sum(meanrep),
               mean = mean(meanrep, na.rm=TRUE),
               median=median(meanrep),
               max=max(meanrep),
               min= min(meanrep),
               sd   = sd(meanrep, na.rm=TRUE),
               se   = sd / sqrt(N)
)
pcrinf$mean<- formatC(pcrinf$mean,  format = "e", digits = 2)
pcrinf$sd<-formatC( pcrinf$sd,  format = "e", digits = 2)
pcrinf$se<-formatC( pcrinf$se,  format = "e", digits = 2)
pcrinf$median<- formatC(pcrtreat$median,  format = "e", digits = 2)
pcrinf$max<-formatC( pcrtreat$max,  format = "e", digits = 2)
pcrinf$min<-formatC( pcrtreat$min,  format = "e", digits = 2)

pcrsum<- ddply(pcr, c("treatment", "day"), summarise,
               N    = sum(meanrep),
               mean = mean(meanrep, na.rm=TRUE),
               median=median(meanrep),
               max=max(meanrep),
               min= min(meanrep),
               sd   = sd(meanrep, na.rm=TRUE),
               se   = sd / sqrt(N)
)
pcrsum$mean<- formatC(pcrsum$mean,  format = "e", digits = 2)
pcrsum$sd<-formatC( pcrsum$sd,  format = "e", digits = 2)
pcrsum$se<-formatC( pcrsum$se,  format = "e", digits = 2)
pcrsum$median<- formatC(pcrtreat$median,  format = "e", digits = 2)
pcrsum$max<-formatC( pcrtreat$max,  format = "e", digits = 2)
pcrsum$min<-formatC( pcrtreat$min,  format = "e", digits = 2)
## mean and median values are very different, hintning to a skewed distribution

### Data visualization

plot1<-ggplot(pcr, aes(x=treatment, y=meanrep, fill=treatment ))+geom_boxplot()+
  ylab(expression("Number of gene copies µL"^-1 ))+
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  scale_fill_manual(values = c("red", "blue", "darkgreen"))+
  theme(text = element_text(size = 16, colour = "black"))+
  theme( axis.text = element_text(colour = "black", size = 16), 
         axis.line = element_line(color="black", size = 1))+
  theme(panel.grid = element_blank(), legend.title = element_blank(), 
        legend.position = "right")+
  facet_grid(~Period)


plot2<-ggplot(pcr, aes(x=Period, y=meanrep, fill=Period ))+ geom_boxplot()+
  ylab(expression("Number of gene copies µL"^-1 ))+ 
  scale_color_manual(values = c("red", "blue", "darkgreen")) + 
  theme(text = element_text(size = 16, colour = "black"))+
  theme( text=element_text(family="Arial", size=16),
         axis.text = element_text(colour = "black", size = 16),
         axis.title.y = element_text(vjust = 2),
         axis.title.x = element_blank(),
         axis.line = element_line(color="black", size = 1),
         panel.grid = element_blank(), legend.title = element_blank(),
         legend.position = "right")+
  facet_grid(~idinf)
plot2
library("gridExtra")
plotpcr<-ggarrange(plot1, plot2, ncol=1, nrow=2)
library(ggpubr)

plotpcr1<-as_ggplot(plotpcr)
highqualgraphR(plot1, "boxplot PCR per treatment", extension = "tiff",
               family="Arial")
highqualgraphR(plot2, "boxplot PCR per day", extension = "tiff",
               family="Arial")
highqualgraphR(plotpcr, "boxplot PCR ", extension = "tiff",
               family="Arial")

#the data are very skewed and therefore are not suitable for tests assuming normal distribution. 
#in addition, sample size is small


##### Make diferent datasets per day and of days we need to compare
day0<-filter(pcr, day==0) ## start experiment and 
day14<-filter(pcr, day==14)## end of stabilization
day15<-filter(pcr, day==15)## first day of treatment
day21<-filter(pcr, day==21)## seventh day of treatment
day27<-filter(pcr, day==27) ## last day of treatment

day0_14<-filter(pcr, day==0 | day==14) ## start experiment and end of stabilization
day0_15<-filter(pcr, day==0| day==15)  ## start experiment and first day of treatment
day14_15<-filter(pcr, day==14| day==15)## end of stabilization and first day of treatment
day14_21<-filter(pcr,  day==14|day==21)## end of stabilization and seventh day of treatment
day14_27<-filter(pcr,  day==14|day==27)## end of stabilization and last day of treatment
pcr1<-filter(pcr, idinfant==1)         ## dataset infant 1
pcr2<-filter(pcr, idinfant==2)         ## dataset infant 2
pcr3<-filter(pcr, idinfant==3)         ## dataset infant 3
datatreat<-filter(pcr, day<15)         ## dataset treatment days



####' as data are skewed with very few observations, we will conduct 
####' non parametric tests: wilcoxon signed rank test and kruskal_Wallis test

## comparing number of gene copies between days 
comp0_14 <- wilcox.test(meanrep ~ day, data = day0_14, paired = TRUE)
comp0_14 #V = 45, p-value = 0.003906
comp14_15<- wilcox.test(meanrep ~ day, data = day14_15, paired = TRUE)
comp14_15 # V = 18, p-value = 0.6523
comp14_21 <- wilcox.test(meanrep ~ day, data = day14_21, paired = TRUE)
comp14_21 #V = 21, p-value = 0.9102
comp14_27 <- wilcox.test(meanrep ~ day, data = day14_27, paired = TRUE)
comp14_27 #V = 37, p-value = 0.09766

## comparing number of gene copies per treatment condition
kruskal.test(meanrep ~ treatment, data = day0)
#Kruskal-Wallis chi-squared = 0.26667, df = 2, p-value = 0.8752
kruskal.test(meanrep ~ treatment, data = day14)
#Kruskal-Wallis chi-squared = 0.8, df = 2, p-value = 0.6703

kruskal.test(meanrep ~ treatment, data = day15)
#Kruskal-Wallis chi-squared = 0.62222, df = 2, p-value = 0.7326
kruskal.test(meanrep ~ treatment, data = day21)
#Kruskal-Wallis chi-squared = 2.7556, df = 2, p-value = 0.2521
kruskal.test(meanrep ~ treatment, data = day27)
#Kruskal-Wallis chi-squared = 1.6889, df = 2, p-value = 0.4298

#################################################################################"
##############################################################################"##
##
##        VFA ANALYSIS
#################################################################################"
##############################################################################"##
##
setwd("C:/Users/latoe/Dropbox/PhD/Labmet/SHIME LNS/1.SHIME results")
vfa<-read.csv('vfa.csv', header = TRUE, sep = ",")

library("ggplot2")
theme_set(theme_bw())
library(extrafont)
#font_import()
#y
loadfonts(device="win") 
library("dplyr")
library("tidyr")
library("ggpubr")
source("C:/Users/latoe/Dropbox/PhD/Labmet/SHIME LNS/1.SHIME results/tryphyloseq/try2/results_whole_treatment_period/Highqualgraphr.r")

#CHECK VARIABLE TYPES
class(vfa$acetate)
class(vfa$propionate)
class(vfa$butyrate)
class(vfa$totalBCFA)
class(vfa$day)
class(vfa$time)
class(vfa$treatment)
class(vfa$period)
class(vfa$idinfant)
vfa$idinf<-factor(vfa$idinfant, labels=c("Infant 1", "Infant 2", "Infant 3"))
vfa$treat<-factor(vfa$treatnum, labels=c("Control", "LNS", "LNSp"))
vfa$treatment<-factor(vfa$treatment, labels=c("Control", "LNS", "LNSp"))
vfa$timef<-factor(vfa$time)
vfa$id<-factor(vfa$idinfant)



dplyr::tbl_df(vfa)


###SUMMARY STATISTICS###

#remove data of time o and27 (no treatment periods)
vfa<-filter(vfa, time>0)
vfa<-filter(vfa, time<14)

#Summary statistics of produced SCFAs in total, then by treatment condition
##summary acetate
sumacet<-vfa %>%
  summarise(mean = mean(acetate),
            sd=sd(acetate),
            median=median(acetate),
            IQR= IQR(acetate),
            min=min(acetate),
            max=max(acetate))
sumacettreat<-vfa %>%
  group_by(treatment) %>%
  summarise(mean = mean(acetate),
            count(vfa),
            sd=sd(acetate),
            median=median(acetate),
            IQR= IQR(acetate),
            min=min(acetate),
            max=max(acetate))
##summary propionate
sumprop<-vfa %>%
  summarise(mean = mean(propionate),
            sd=sd(propionate),
            median=median(propionate),
            IQR= IQR(propionate),
            min=min(propionate),
            max=max(propionate))
sumproptreat<-vfa %>%
  group_by(treatment) %>%
  summarise(mean = mean(propionate),
            sd=sd(propionate),
            median=median(propionate),
            IQR= IQR(propionate),
            min=min(propionate),
            max=max(propionate))


##summary butyrate
sumbut<-vfa %>%
  summarise(mean = mean(butyrate),
            sd=sd(butyrate),
            median=median(butyrate),
            IQR= IQR(butyrate),
            min=min(butyrate),
            max=max(butyrate))
sumbuttreat<-vfa %>%
  group_by(treatment) %>%
  summarise(mean = mean(butyrate),
            sd=sd(butyrate),
            median=median(butyrate),
            IQR= IQR(butyrate),
            min=min(butyrate),
            max=max(butyrate))
##summary BCFA
sumbc<-vfa %>%
  summarise(mean = mean(totalBCFA),
            sd=sd(totalBCFA),
            median=median(totalBCFA),
            IQR= IQR(totalBCFA),
            min=min(totalBCFA),
            max=max(totalBCFA))
sumbctreatt<-vfa %>%
  group_by(treatment) %>%
  summarise(mean = mean(totalBCFA),
            sd=sd(totalBCFA),
            median=median(totalBCFA),
            IQR= IQR(totalBCFA),
            min=min(totalBCFA),
            max=max(totalBCFA))

#Checking the data distribution
##histograms
#graphics.off()
par("mar")
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
hist(vfa$acetate, breaks = 28, main = "Acetate", xlab = NULL)
hist(vfa$propionate, breaks = 28, main = "Propionate", xlab = NULL)
hist(vfa$butyrate, breaks = 28, main = "Butyrate", xlab = NULL)
hist(vfa$totalBCFA, breaks = 28, main = "BCFA", xlab = NULL)
## Density plots
dens1<-density(vfa$acetate)
dens2<-density(vfa$propionate)
dens3<-density(vfa$butyrate)
dens4<-density(vfa$totalBCFA)

#graphics.off()
par("mar")
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
plot(dens1, frame=FALSE, col="blue", main = "acetate", cex.main=0.75, cex.axis=0.75, cex.lab=0.75)
plot(dens2, frame=FALSE, col="blue", main = "propionate ", cex.main=0.75, cex.axis=0.75, cex.lab=0.75)
plot(dens3, frame=FALSE, col="blue", main = "butyrate" , cex.main=0.75, cex.axis=0.75, cex.lab=0.75)
plot(dens4, frame=FALSE, col="blue", main = "BCFA", cex.main=0.75, cex.axis=0.75, cex.lab =0.75)
### Data are skewed and slightly bimodal for acetate, propionate an butyrate. 
### BCFA is symetrically distributed

#### Scatterplots of SCFA production per infant 
par("mar")
par(mar=c(4,4,4,4))
scatacet<-ggplot(vfa, aes(x=day, y=acetate, col=treatment))+geom_point(size=3)+ geom_line()+
  xlab("Day")+ ylab(expression("Acetate mmol L"^-1))+
  scale_color_manual(name= "Treatment",breaks=c("Control","LNS", "LNSp"),values = c("red", "blue", "darkgreen")) +
  scale_x_continuous(breaks = round(seq(min(vfa$day), max(vfa$day), by = 3),1))+
  theme(text=element_text(family="Arial", size=16))+
  theme(axis.title.y = element_text(vjust = 2), axis.text = element_text(colour = "black", size = 16))+
  facet_grid(~idinf)+theme(panel.grid = element_blank())


scatacet
scatprop<-ggplot(vfa, aes(x=day, y=propionate, col=treatment))+geom_point(size=3)+ geom_line()+
  xlab("Day")+ ylab(expression("Propionate mmol L"^-1))+ 
  scale_x_continuous(breaks = round(seq(min(vfa$day), max(vfa$day), by = 3),1))+
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  theme(text=element_text(family="Arial", size=16))+
  theme(axis.title.y = element_text(vjust = 2), axis.text = element_text(colour = "black", size = 16))+
  facet_grid(~idinf)+theme(panel.grid = element_blank())

scatbut<-ggplot(vfa, aes(x=day, y=butyrate, col=treatment))+geom_point(size=3)+ geom_line()+
  xlab("Day")+ ylab(expression("Butyrate mmol L"^-1))+
  scale_x_continuous(breaks = round(seq(min(vfa$day), max(vfa$day), by = 3),1))+
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  theme(text=element_text(family="Arial", size=16))+
  theme(axis.title.y = element_text(vjust = 2), axis.text = element_text(colour = "black", size = 16))+
  facet_grid(~idinf)+theme(panel.grid = element_blank())
scatbut
scatbc<-ggplot(vfa, aes(x=day, y=totalBCFA, col=treatment))+geom_point(size=3)+ geom_line()+
  xlab("Day")+ ylab(expression("BCFA mmol L"^-1))+ 
  scale_x_continuous(breaks = round(seq(min(vfa$day), max(vfa$day), by = 3),1))+
  scale_color_manual("Treatment",values = c("red", "blue", "darkgreen")) +
  theme(text=element_text(family="Arial", size=16))+
  theme(axis.title.y = element_text(vjust = 2), axis.text = element_text(colour = "black", size = 16))+
  facet_grid(~idinf)+theme(panel.grid = element_blank())


library("gridExtra")
scatplotvfa<-grid.arrange(scatacet, scatprop, scatbut, scatbc,  ncol=2, nrow=2)
p5<-as_ggplot(scatplotvfa)
highqualgraphR(p5,"figure 5",extension = "tiff",
               family="Arial")

######??? correlation between OTU and SCFA
## get the OTU tables in vegan
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTU1<-veganotu(shime1)
OTU2<-veganotu(shime2)
OTU3<-veganotu(shime3)
## format environmental data per infant
env.data <- filter(map0, Day>14 & Day<28)
env.data$Day<-factor(env.data$Day)
#labels = c("Day15", "Day16", "Day17",
#          "Day18", "Day19", "Day20",
#         "Day21", "Day22", "Day23",
#        "Day24", "Day25", "Day26", "Day27"
## create a color column for treatment
env.data$col="red"
env.data$col[env.data$Treatment=="LNS"]="blue"
env.data$col[env.data$Treatment=="LNSp"]="darkgreen"
rownames(env.data)=env.data$SampleName
env.data = env.data[,-which(names(env.data) %in% c("SampleID", "SampleName"))]

#make environmental data for each infant
env.data1 <- subset(env.data, idinfant=="Infant 1")
env.data2 <- subset(env.data, idinfant=="Infant 2")
env.data3 <- subset(env.data, idinfant=="Infant 3")
#??? extract the taxonomy data
tax1<-as.data.frame(shime1@tax_table)
tax2<-as.data.frame(shime2@tax_table)
tax3<-as.data.frame(shime3@tax_table)
#' assign genus as column names  and sampleNames 
#' as rownames in OTU tables
colnames(OTU1)=tax1$Genus
rownames(OTU1)=rownames(env.data1)
colnames(OTU2)=tax2$Genus
rownames(OTU2)=rownames(env.data2)
colnames(OTU3)=tax3$Genus
rownames(OTU3)=rownames(env.data3)

#Calculate relative abundances to be used to select the most relevant OTUs

abund1 = OTU1/rowSums(OTU1)*100
abund2 = OTU2/rowSums(OTU2)*100
abund3 = OTU3/rowSums(OTU3)*100

#Subset OTUs to abundance cutoff (OTUs with more that 1% abundance are selected)
OTU.abund1 = OTU1[, apply(abund1, MARGIN=2, function(x) any(x > 1))]
OTU.abund2 = OTU2[, apply(abund2, MARGIN=2, function(x) any(x > 1))]
OTU.abund3 = OTU3[, apply(abund3, MARGIN=2, function(x) any(x > 1))]

cor.pearson1ac = as.data.frame(cor(OTU.abund1, env.data1$acetate, method = "pearson"))
cor.pearson1pr = as.data.frame(cor(OTU.abund1, env.data1$propionate, method = "pearson"))
cor.pearson1bu = as.data.frame(cor(OTU.abund1, env.data1$butyrate, method = "pearson"))
cor.pearson1bc = as.data.frame(cor(OTU.abund1, env.data1$BCFA, method = "pearson"))

names(cor.pearson1ac)[1] <- "Acetate"
names(cor.pearson1pr)[1] <- "Propionate"
names(cor.pearson1bu)[1] <- "Butyrate"
names(cor.pearson1bc)[1] <- "BCFA"
ac1<-subset(cor.pearson1ac, Acetate>0.5)
pr1<-subset(cor.pearson1pr, Propionate>0.5)
bu1<-subset(cor.pearson1bu, Butyrate>0.5)
bc1<-subset(cor.pearson1bc, BCFA>0.5)


cor.pearson2ac = as.data.frame(cor(OTU.abund2, env.data1$acetate, method = "pearson"))
cor.pearson2pr = as.data.frame(cor(OTU.abund2, env.data1$propionate, method = "pearson"))
cor.pearson2bu = as.data.frame(cor(OTU.abund2, env.data1$butyrate, method = "pearson"))
cor.pearson2bc = as.data.frame(cor(OTU.abund2, env.data1$BCFA, method = "pearson"))

names(cor.pearson2ac)[1] <- "Acetate"
names(cor.pearson2pr)[1] <- "Propionate"
names(cor.pearson2bu)[1] <- "Butyrate"
names(cor.pearson2bc)[1] <- "BCFA"
ac2<-subset(cor.pearson2ac, Acetate>0.5)
pr2<-subset(cor.pearson2pr, Propionate>0.5)
bu2<-subset(cor.pearson2bu, Butyrate>0.5)
bc2<-subset(cor.pearson2bc, BCFA>0.5)

cor.pearson3ac = as.data.frame(cor(OTU.abund3, env.data1$acetate, method = "pearson"))
cor.pearson3pr = as.data.frame(cor(OTU.abund3, env.data1$propionate, method = "pearson"))
cor.pearson3bu = as.data.frame(cor(OTU.abund3, env.data1$butyrate, method = "pearson"))
cor.pearson3bc = as.data.frame(cor(OTU.abund3, env.data1$BCFA, method = "pearson"))

names(cor.pearson3ac)[1] <- "Acetate"
names(cor.pearson3pr)[1] <- "Propionate"
names(cor.pearson3bu)[1] <- "Butyrate"
names(cor.pearson3bc)[1] <- "BCFA"
ac3<-subset(cor.pearson3ac, Acetate>0.5)
pr3<-subset(cor.pearson3pr, Propionate>0.5)
bu3<-subset(cor.pearson3bu, Butyrate>0.5)
bc3<-subset(cor.pearson3bc, BCFA>0.5)

# INFERENTIAL ANALYSIS
## CHECKING ASSUMPTIONS FOR ANALYSIS
### check normality of observations  with QQ-plots
graphics.off()
par("mar")
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
qqnorm(vfa$acetate, pch=1,main = "Normal Q-Qplot acetate", cex.lab=0.75, cex.main=0.75)
qqline(vfa$acetate, col = "red", lwd = 2)

qqnorm(vfa$propionate, pch=1, main = "Normal Q-Qplot propionate", cex.lab=0.75, cex.main=0.75)
qqline(vfa$propionate, col = "red", lwd = 2)

qqnorm(vfa$butyrate, pch=1, main = "Normal Q-Qplot butyrate", cex.lab=0.75, cex.main=0.75)
qqline(vfa$butyrate, col = "red", lwd = 2)

qqnorm(vfa$totalBCFA, pch=1,main = "Normal Q-Qplot BCFA", cex.lab=0.75, cex.main=0.75)
qqline(vfa$totalBCFA, col = "red", lwd = 2)
#### Acetate, propionate and butyrate deviate from normal distribution, BCFA seems normally distributed 

###Shapiro wilk test of normality 
shapiro.test(vfa$acetate)
shapiro.test(vfa$propionate)
shapiro.test(vfa$butyrate)
shapiro.test(vfa$totalBCFA)
#### the Shapiro-Wilk's test confirms what was oberved with QQ-plots

### Repeat shapiro wilk test per treatment condition
tests = vfa %>%
  group_by(treatment) %>%
  do(test = shapiro.test(.$acetate))
tests$test

vfa %>%
  group_by(treatment) %>%
  summarise(stest = shapiro.test(acetate)$p.value)

vfa %>%
  group_by(treatment) %>%
  summarise(stest = shapiro.test(propionate)$p.value)

vfa %>%
  group_by(treatment) %>%
  summarise(stest = shapiro.test(butyrate)$p.value)

vfa %>%
  group_by(treatment) %>%
  summarise(stest = shapiro.test(totalBCFA)$p.value)
#### all SCFA production follow normal distribution within treatment conditions
#### except butyrate production in LNSp treatment condition

### Repeat shapiro wilk test per infant
vfa %>%
  group_by(idinf) %>%
  summarise(stest = shapiro.test(acetate)$p.value)

vfa %>%
  group_by(idinf) %>%
  summarise(stest = shapiro.test(propionate)$p.value)

vfa %>%
  group_by(idinf) %>%
  summarise(stest = shapiro.test(butyrate)$p.value)

vfa %>%
  group_by(idinf) %>%
  summarise(stest = shapiro.test(totalBCFA)$p.value)
#### Actate, propionate and butyrate production are not normally distruted within infants. 
#### BCFA production deviates from normal distribution in infant 1

### Repeat shapiro wilk test per infant and treatment condition
vfa %>%
  filter(idinf== "Infant 1")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(acetate)$p.value)
vfa %>%
  filter(idinf== "Infant 2" )%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(acetate)$p.value)

vfa %>%
  filter(idinf== "Infant 3")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(acetate)$p.value)

vfa %>%
  filter(idinf== "Infant 1")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(propionate)$p.value)
vfa %>%
  filter(idinf== "Infant 2" )%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(propionate)$p.value)

vfa %>%
  filter(idinf== "Infant 3")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(propionate)$p.value)

vfa %>%
  filter(idinf== "Infant 1")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(butyrate)$p.value)
vfa %>%
  filter(idinf== "Infant 2" )%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(butyrate)$p.value)

vfa %>%
  filter(idinf== "Infant 3")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(butyrate)$p.value)

vfa %>%
  filter(idinf== "Infant 1")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(totalBCFA)$p.value)
vfa %>%
  filter(idinf== "Infant 2" )%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(totalBCFA)$p.value)

vfa %>%
  filter(idinf== "Infant 3")%>%
  group_by(treatment)%>%
  summarise(stest = shapiro.test(totalBCFA)$p.value)
#### not normally distributed:  
#### butyrate in LNSp, infant 2 and 3
#### propionate in LNSp, infant2


#check linear relationship with Time, per treatment group and visualizing the slopes and intercepts per treatment and infant
##Acetate
gacet<-ggscatter(vfa, x = "day", y = "acetate", color = c("treatment"),
                 palette=c("red", "blue", "darkgreen"),
                 add = "reg.line", conf.int = TRUE,
                 xlab = "Sampling tday", ylab = "Acetate (mMol/L)")
gacet<-gacet+scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  theme(
    legend.position = "bottom",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())
##propionate
gprop<-ggscatter(vfa, x = "day", y = "propionate", color = "treatment",
                 palette=c("red", "blue", "darkgreen"),
                 add = "reg.line", conf.int = TRUE, 
                 xlab = "Sampling day", ylab = "Propionate (mMol/L)")
gprop<-gprop+scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  theme(
    legend.position = "bottom",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())

##butyrate
gbut<-ggscatter(vfa, x = "day", y = "butyrate", color = "treatment",
                palette=c("red", "blue", "darkgreen"),
                add = "reg.line", conf.int = TRUE, 
                xlab = "Sampling day", ylab = "Butyrate (mMol/L)")
gbut<-gbut+scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  theme(
    legend.position = "bottom",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())

##BCFA
gbc<-ggscatter(vfa, x = "day", y = "totalBCFA", color = "treatment",
               palette=c("red", "blue", "darkgreen"),
               add = "reg.line", conf.int = TRUE, 
               xlab = "Sampling day", ylab = "BCFA (mMol/L)")
gbc<-gbc+scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  theme(text=element_text(family="Arial", size=12),
        legend.position = "bottom",
        axis.title= element_text(color="black", size=12, face="bold"),
        legend.title = element_blank())
scattervfa<-ggarrange(gacet, gprop, gbut, gbc,  ncol=2, nrow=2, label.x=0,common.legend = TRUE, legend = "right")
annotate_figure(scattervfa, bottom = text_grob("Short Chain Fatty Acids production over day per infant"))
###' the relation between SCFA production and day is not linear. slopes of production do not show a common pattern per 
###' treatment group but seem to depend on infant. intercepts seem to be different per treatment and infant

#visualization with scatter plot and Loess smoothers
##Acetate
gacet1<-ggscatter(vfa, x = "day", y = "acetate", color = "treatment",
                  palette=c("red", "blue", "darkgreen"),
                  add = "loess", conf.int = TRUE,
                  xlab = "Sampling day", ylab = "Acetate (mMol/L)")

gacet1<-gacet1+
  scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  ggtitle("Acetate production")+
  theme(
    legend.position = "none",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())
gacet1
##propionate
gprop1<-ggscatter(vfa, x = "day", y = "propionate", color = "treatment",
                  palette=c("red", "blue", "darkgreen"),
                  add = "loess", conf.int = TRUE,
                  xlab = "Sampling day", ylab = "Propionate (mMol/L)")
gprop1<-gprop1+
  scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  ggtitle("Propionate production" )+
  theme(
    legend.position = "none",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())
gprop1
##butyrate
gbut1<-ggscatter(vfa, x = "day", y = "butyrate", color = "treatment",
                 palette=c("red", "blue", "darkgreen"),
                 add = "loess", conf.int = TRUE,
                 xlab = "Sampling day", ylab = "Butyrate (mMol/L)")
gbut1<-gbut1+
  scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  ggtitle("Butyrate production")+
  theme(
    legend.position = "none",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())
gbut1
##BCFA
gbc1<-ggscatter(vfa, x = "day", y = "totalBCFA", color = "treatment",
                palette=c("red", "blue", "darkgreen"),
                add = "loess", conf.int = TRUE,
                xlab = "Sampling day", ylab = "BCFA (mMol/L)")
gbc1<-gbc1+
  scale_x_continuous(breaks = pretty(vfa$day, n = 5)) + facet_grid(.~idinf)+
  ggtitle("BCFA production")+
  theme(
    legend.position = "none",
    plot.title = element_text(color="black", size=12, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_blank())

scattervfa1<-ggarrange(gacet1, gprop1, gbut1, gbc1,  ncol=2, nrow=2, label.x=0,common.legend = TRUE, legend = "right")
annotate_figure(scattervfa1, bottom = text_grob("Short Chain Fatty Acids production over day per infant"))
### Loess smoothers describe better the relation between SCFA production and day
### The intercepts seem to depend on the combination of infant and treatment, and not on infant or treatment alone

############################
############################  Using lmer to define models for inferential analysis  ####

library("lme4")

library("emmeans")
library("lmerTest")
library("lattice")


####BUILDING MODELS##########
#'fitting full models with Vector-valued random effects (random slope of treatment and random intercept of infant) by subject and 
#'reduced models with Random intercepts for infant  and infant /treatment (scalar)
#'this is based on the visualization
#'for eah outcome, model 2 and 3 have the same number of parameters and are not nested can only be compared using AIC

##  Model for acetate production  ##
###################################

ma0<-lmer(acetate~treat+day+(treat|idinf),REML=FALSE, data=vfa) #singular fit. the model may be overfitted and will not be further considered
ma1<-lmer(acetate~treat+day +(1|idinf)+(1|idinf:treat),REML=FALSE, data=vfa) 
ma2<-lmer(acetate~treat+day +(1|idinf:treat),REML=FALSE, data=vfa)
ma3<-lmer(acetate~treat+day +(1|idinf),REML=FALSE, data=vfa)
AIC(ma0, ma1, ma2, ma3)
anova(ma0, ma1, ma2, ma3)
(pb <- pbkrtest::PBmodcomp(ma0,ma1,seed=101)) 
#Parametric bootstrap test; time: 100.90 sec; samples: 1000 extremes: 300;
#Requested samples: 1000 Used samples: 995 Extremes: 300
#large : acetate ~ treat + day + (treat | idinf)
#small : acetate ~ treat + day + (1 | idinf) + (1 | treat:idinf)
#stat df p.value
#LRT    4.0361  4  0.4011
#PBtest 4.0361     0.3022
(pb <- pbkrtest::PBmodcomp(ma1,ma2,seed=101)) 
#Parametric bootstrap test; time: 52.38 sec; samples: 1000 extremes: 15;
#Requested samples: 1000 Used samples: 608 Extremes: 15
#large : acetate ~ treat + day + (1 | idinf) + (1 | treat:idinf)
#small : acetate ~ treat + day + (1 | treat:idinf)
#stat df  p.value   
#LRT    7.4507  1 0.006341 **
#PBtest 7.4507    0.026273 * 
(pb <- pbkrtest::PBmodcomp(ma2,ma3,seed=101)) 

#'Parametric bootstrap test; time: 41.83 sec; samples: 1000 extremes: 29;
#'Requested samples: 1000 Used samples: 29 Extremes: 29
#'large : acetate ~ treat + day + (1 | treat:idinf)
#'small : acetate ~ treat + day + (1 | idinf)
#'           stat df p.value
#'LRT    -6.5203  0       1
#'PBtest -6.5203          1
(pb <- pbkrtest::PBmodcomp(ma1,ma3,seed=101))
#Parametric bootstrap test; time: 50.80 sec; samples: 1000 extremes: 46;
#Requested samples: 1000 Used samples: 478 Extremes: 46
#large : acetate ~ treat + day + (1 | idinf) + (1 | treat:idinf)
#small : acetate ~ treat + day + (1 | idinf)
#stat df p.value  
#LRT    0.9304  1 0.33475  
#PBtest 0.9304    0.09812 .
anova (ma1, ma3)
### THE FINAL MODEL IS  MA3 based on AIC


##   Model for propionate production  ##
#######################################

mp0<-lmer(propionate~treat+day+(treat|idinf),REML=FALSE, data=vfa) #singular fit. the model may be  overfitted 
#and will not be further considered
mp1<-lmer(propionate~treat+day +(1|idinf)+(1|idinf:treat),REML=FALSE, data=vfa) 
mp2<-lmer(propionate~treat+day +(1|idinf:treat),REML=FALSE, data=vfa)
mp3<-lmer(propionate~treat+day +(1|idinf),REML=FALSE, data=vfa)
AIC(mp0, mp1, mp2, mp3)
anova(mp0, mp1, mp2, mp3)

(pb <- pbkrtest::PBmodcomp(mp1,mp2,seed=101)) 
#Parametric bootstrap test; time: 60.04 sec; samples: 1000 extremes: 110;
#Requested samples: 1000 Used samples: 634 Extremes: 110
#large : propionate ~ treat + day + (1 | idinf) + (1 | idinf:treat)
#small : propionate ~ treat + day + (1 | idinf:treat)
#stat df p.value
#LRT    1.8087  1  0.1787
#PBtest 1.8087     0.1748

(pb <- pbkrtest::PBmodcomp(mp2,mp3,seed=101)) 
#Parametric bootstrap test; time: 41.86 sec; samples: 1000 extremes: 0;
#Requested samples: 1000 Used samples: 37 Extremes: 0
#large : propionate ~ treat + day + (1 | idinf:treat)
#small : propionate ~ treat + day + (1 | idinf)
#       stat    df   p.value    
#LRT    15.752  0    < 2e-16 ***
#PBtest 15.752       0.02632 *  
(pb <- pbkrtest::PBmodcomp(mp1,mp3,seed=101)) 
#Parametric bootstrap test; time: 54.52 sec; samples: 1000 extremes: 0;
#Requested samples: 1000 Used samples: 459 Extremes: 0
#large : propionate ~ treat + day + (1 | idinf) + (1 | idinf:treat)
#small : propionate ~ treat + day + (1 | idinf)
#         stat df   p.value    
#LRT    17.561  1   2.782e-05 ***
#PBtest 17.561     0.002174 ** 

### THE FINAL MODEL IS MP2

## Model for butyrate production ##
##################################

mb0<-lmer(butyrate~treat+day+(treat|idinf),REML=FALSE, data=vfa) #singular fit and does not coverge. the model may be  overfitted
mb1<-lmer(butyrate~treat+day +(1|idinf)+(1|idinf:treat),REML=FALSE, data=vfa) # singular fit
mb2<-lmer(butyrate~treat+day +(1|idinf:treat),REML=FALSE, data=vfa)
mb3<-lmer(butyrate~treat+day +(1|idinf),REML=FALSE, data=vfa)# singular fit
## Only model 2 seems to work. try a bayesian method to avoid singularity 
#install.packages("blme")
library(blme)
mb00<-blmer(butyrate~treat +time+(treat|idinf),REML=FALSE, data=vfa)
mb01<-blmer(butyrate~treat +time+(1|idinf)+(1|treat:idinf),REML=FALSE, data=vfa) #" does not converge
mb02<-blmer(butyrate~treat +time+(1|idinf),REML=FALSE, data=vfa)
mb03<-blmer(butyrate~treat +time+(1|idinf:treat),REML=FALSE, data=vfa) # does not converge

AIC(mb00, mb01, mb02, mb03)
anova(mb00, mb01, mb02, mb03)
(pb <- pbkrtest::PBmodcomp(mb00,mb02,seed=101)) #NOT APPLICABLE DF=0
# MB2 IS NOT OVERFIT AND HAS THE BEST AIC
#the final model is thus mb2

##  Model for BCFA production  ##
################################

mbc0<-lmer(totalBCFA~treat+day+(treat|idinf),REML=FALSE, data=vfa) #singular fit and no convergence. the model may be  overfitted and will not be further considered
mbc1<-lmer(totalBCFA~treat+day +(1|idinf)+(1|idinf:treat),REML=FALSE, data=vfa) # singular fit 
mbc2<-lmer(totalBCFA~treat+day +(1|idinf:treat),REML=FALSE, data=vfa)
mbc3<-lmer(totalBCFA~treat+day +(1|idinf),REML=FALSE, data=vfa)

(pb <- pbkrtest::PBmodcomp(mbc2,mbc3,seed=101)) #not applicabel df=0 (same number of terms)


AIC(mbc2, mbc3)
anova(mbc2, mbc3)
# THE FINAL MODEL IS mbc3 BASED ON AIC


###Final Models
##    final model for acetate production

acet<-update(ma3, REML=TRUE)
prop<-update(mp2, REML=TRUE)
but<-update(mb2, REML=TRUE)
bcfa<-update(mbc3, REML=TRUE)

### post-hoc check of models                                                                                                           

#graphics.off()
par("mar")
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
plot(acet, main="Acetate")
plot(prop, main="Propionate")
plot(but, main="Butyate")
plot(bcfa, main="BCFA")

graphics.off()
par("mar")
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))

qqnorm(resid(acet), main = "Residual QQplot acetate", cex=2, pch=16);qqline(resid(acet), lwd=4)

qqnorm(resid(prop), main = "Residual QQplot propionate", cex=2, pch=16)
qqline(resid(prop), lwd=4)

qqnorm(resid(but), main = "Residual QQplot butyrate", cex=2, pch=16)
qqline(resid(but), lwd=4)

qqnorm(resid(bcfa), main = "Residual QQplot BCFA", cex=2, pch=16)
qqline(resid(bcfa), lwd=4)
#### post hoc visualization shows that residuals of all model do not severely deviate from normality,
#### thus validating the models

#### Models summary and contrasts using Tukey correction for multiple comparisons
summary(acet)
summary(prop)
summary(but)
summary(bcfa)
######  contrasts xith tukey correction
pacet<-emmeans(acet, "treat", adjust="Tukey")
pacet
pprop<-emmeans(prop, "treat", adjust="Tukey")
pprop
pbut<-emmeans(but,"treat", adjust="Tukey")
pbcfa<-emmeans(bcfa, "treat", adjust="Tukey")
pbut
pbcfa
plot(pacet, comparisons = TRUE)
plot(pprop, comparisons = TRUE)
plot(pbut, comparisons = TRUE)
plot(pbcfa, comparisons = TRUE)

pmmeans(acet, revpairwise~treat, adjust="Tukey")
pmmeans(prop, revpairwise~treat, adjust="Tukey")
pmmeans(but, revpairwise~treat, adjust="Tukey")
pmmeans(bcfa, revpairwise~treat, adjust="Tukey")


