### This code analyzes the pilot study data
### 1. DNA concentrations in different sample preparations
### 2. DNA sequence analysis: read counts & diversity in different sample preparations
### Last updated 5.23.2023 by Kara Andres (akara@wustl.edu)

### Clear work environment and load packages
rm(list = ls())
library(ggplot2)
library(phyloseq)
library(qiime2R)
library(reshape2)
library(decontam)
library(Biostrings)
library(stringr)
library(phylosmith)
library(RColorBrewer)
library(colorspace)
library(forcats)
library(vegan)
library(ggpubr)
library(rstatix)
library(pander)
library(dplyr)
library(lefser)

##############################################################################
############ DNA concentrations in different sample preparations #############
##############################################################################

### Load dataset
pilot_study_metadata <- read.csv("pilot_study_metadata.csv", header=TRUE)
pilot_study_metadata <- pilot_study_metadata[!pilot_study_metadata$Individual=="blank",]
pilot_study_metadata$Fresh_frozen <- factor(pilot_study_metadata$Fresh_frozen, levels=c("fresh", "frozen"), labels=c("Fresh", "Frozen"))
pilot_study_metadata$Extraction <- factor(pilot_study_metadata$Extraction, levels=c("washed", "direct"), labels=c("Washed", "Unwashed"))
pilot_study_metadata$Species <- factor(pilot_study_metadata$Species, levels=c("Campostoma anomalum", "Etheostoma caeruleum"), labels=c("Herbivore", "Insectivore"))
species.cols <- c("#82C341", "#FDB81A")
fh.fz.cols <- c("salmon", "lightblue")

### Grouped histogram of fish lengths
ggplot(pilot_study_metadata, aes(x=Length_cm, fill=Species)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=species.cols) + 
  theme_bw()

### Grouped boxplot of DNA concentrations
conc_plot <- ggplot(pilot_study_metadata, aes(x=Fresh_frozen, y=DNA_concent_ng_ul, fill=Fresh_frozen)) + 
  geom_boxplot() + 
  facet_wrap(~Species + Extraction, scales="free_x") +
  scale_fill_manual(values=fh.fz.cols) +
  labs(x="", y="DNA Concentration (ng/µL)", size=12) + theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text.x = element_text(size=12),
        legend.text = element_text(size=12))
conc_plot
# ggsave("Figures/DNA_concentration.pdf", plot=conc_plot)

### Grouped boxplot of Cq
cq_plot <- ggplot(pilot_study_metadata, aes(x=Fresh_frozen, y=Cq, fill=Fresh_frozen)) + 
  geom_boxplot() + 
  facet_wrap(~Species + Extraction, scales="free_x") +
  scale_fill_manual(values=fh.fz.cols) +
  labs(x="", y="Cq", size=12) + theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text.x = element_text(size=12),
        legend.text = element_text(size=12))
cq_plot
# ggsave("Figures/Cq_values.pdf", plot=cq_plot)

aggregate(x=pilot_study_metadata$Cq, by=list(pilot_study_metadata$Fresh_frozen), FUN=mean, na.rm=TRUE)
aggregate(x=pilot_study_metadata$Cq, by=list(pilot_study_metadata$Extraction), FUN=mean, na.rm=TRUE)
aggregate(x=pilot_study_metadata$Cq, by=list(pilot_study_metadata$Species), FUN=mean, na.rm=TRUE)

##############################################################################
############### Data pruning, decontamination, and rarefaction ###############
##############################################################################

### Produce phyloseq object
SVs <- read_qza("table.qza")
taxonomy <- read_qza("taxonomy.qza")
metadata <- read_q2metadata("meta_data_R.txt")
rownames(metadata) <- metadata$SampleID
metadata$fish_S_ID <- as.factor(metadata$fish_S_ID)
metadata$fh_fz <- as.factor(metadata$fh_fz)
metadata$fish_species <- as.factor(metadata$fish_species)
metadata$washed <- as.factor(metadata$washed)
sampledata <- sample_data(metadata)
rownames(sampledata) <- rownames(metadata)
R_tree <- read_qza("rooted-tree.qza")
phyloseq_dat <- phyloseq(otu_table(SVs$data,taxa_are_rows = TRUE),tax_table(as.matrix(parse_taxonomy(taxonomy$data))),sampledata,phy_tree(R_tree$data))
# save(phyloseq_dat,file = "phyloseq_dat.RDATA")

### Basic denoising 
load("phyloseq_dat.RDATA")
pilot_nontree <- phyloseq_dat
sample_names(pilot_nontree) <- as.character(pilot_nontree@sam_data$SampleID)
df <- data.frame(Sample=sample_names(pilot_nontree), reads=sample_sums(pilot_nontree@otu_table))

### Inspect and remove negative controls
pilot_blanks <- subset_samples(pilot_nontree,fish_species=="blank")
pilot_blanks <- prune_taxa(taxa_sums(pilot_blanks) > 0, pilot_blanks)
SampleFilter1 <- sample_names(pilot_blanks)[sample_sums(pilot_blanks) < 50]
pilot_nontree <- subset_samples(pilot_nontree,!sample_names(pilot_nontree) %in% SampleFilter1)
pilot_nontree_m <- prune_taxa(taxa_sums(pilot_nontree)>3,pilot_nontree)
KingdomFilter <- c("d__Eukaryota","Unassigned")
pilot_nontree_m <- subset_taxa(pilot_nontree_m,!Kingdom %in% KingdomFilter)
# save(pilot_nontree_m,file = "pilot_before_decontam.RDATA")

### Normalize reads in each sample (transform to relative abundance)
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(10), rainbow(10), "azure3")
normld_pilot.ReAb <- transform_sample_counts(pilot_nontree_m, function(OTU) OTU*100/sum(OTU))
base.phylum <- normld_pilot.ReAb %>% tax_glom(taxrank="Phylum", NArm=FALSE) %>% psmelt()

### Calculate correlation between technical replicates
pilot_nontree_m.ReAb <- transform_sample_counts(pilot_nontree_m, function(OTU) OTU/sum(OTU))
Pearson_cor_table <- as.data.frame(cor(pilot_nontree_m@otu_table, use="all.obs", method="pearson"))
keep <- sort(colnames(Pearson_cor_table))
Pearson_cor_table_srt <- Pearson_cor_table[,keep]
Pearson_cor_table_srt <- Pearson_cor_table_srt[order(row.names(Pearson_cor_table_srt)), ]
df <- data.frame(Samplename=NA,Pearson_values=NA)
number_index <- seq(1,96, by=2)
for(i in number_index){
  df[i,2] <- Pearson_cor_table_srt[i, i+1]
  df[i,1] <- rownames(Pearson_cor_table_srt)[i]
}
pearson_table <- na.omit(df)
mean(pearson_table$Pearson_values)

### Identify and remove contaminants by decontam
load("pilot_before_decontam.RDATA")
sample_data(pilot_nontree_m)$is.neg <- sample_data(pilot_nontree_m)$fish_species == "blank"
contamdf.prev05 <- isContaminant(pilot_nontree_m, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
which(contamdf.prev05$contaminant)
decontam_Contaminant <- rownames(contamdf.prev05)[which(contamdf.prev05$contaminant)]
# No contaminants found
# pilot_nontree_m_noncontam <- prune_taxa(!rownames(pilot_nontree_m@otu_table) %in% contaminants,pilot_nontree_m)
# save(pilot_nontree_m,file = "pilot_after_decontam.RDATA")

### Sample rarefaction
load("pilot_after_decontam.RDATA")
sample_data(pilot_nontree_m)$fish_S_ID <- as.factor(sample_data(pilot_nontree_m)$fish_S_ID) # merge replicates
pilot_nontree_m.nonrep <- merge_samples(pilot_nontree_m,"fish_S_ID")

### Metadata coerced to NA -- manually add back sample data
sample_table <- unique(sample_data(pilot_nontree_m)[,c("fish_S_ID","fh_fz","fish_species","washed")])
rownames(sample_table) <- sample_table$fish_S_ID
sample_data(pilot_nontree_m.nonrep) <- sample_table

### Look at sequencing depth by sample type 
read_depth <- data.frame(fish_S_ID=sample_data(pilot_nontree_m.nonrep)$fish_S_ID,
                         fh_fz=sample_data(pilot_nontree_m.nonrep)$fh_fz,
                         fish_species=sample_data(pilot_nontree_m.nonrep)$fish_species,
                         washed=sample_data(pilot_nontree_m.nonrep)$washed,
                         read_depth=sample_sums(pilot_nontree_m.nonrep))
read_depth$fh_fz <- factor(read_depth$fh_fz, levels=c("fresh", "frozen"), labels=c("Fresh", "Frozen"))
read_depth$washed <- factor(read_depth$washed, levels=c("washed", "unwashed"), labels=c("Washed", "Unwashed"))
read_depth$fish_species <- factor(read_depth$fish_species, levels=c("C_anomalum", "E_caeruleum"), labels=c("Herbivore", "Insectivore"))
aggregate(x=read_depth$read_depth, by=list(read_depth$fh_fz), FUN=mean)
aggregate(x=read_depth$read_depth, by=list(read_depth$washed), FUN=mean)
aggregate(x=read_depth$read_depth, by=list(read_depth$fish_species), FUN=mean)
data_hline <- data.frame(aggregate(x=read_depth$read_depth,
          by=list(read_depth$fh_fz,read_depth$washed, read_depth$fish_species),
          FUN=mean))
names(data_hline) <- c("fh_fz", "washed", "fish_species", "hline")
plot1 <- ggplot(data=read_depth, aes(x=fish_S_ID, y=read_depth)) +
  geom_bar(stat="identity") +
  ylab("Read Depth") + xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  geom_hline(data=data_hline, aes(yintercept=hline), color="#F8766D", linetype="dashed") +
  facet_wrap(~fish_species + fh_fz + washed, scales="free_x", nrow=1, ncol=6)
plot1
# ggsave("Figures/read_depth_all_samples.pdf", plot=plot1)

### Remove samples with reads < 300
sample_sums(pilot_nontree_m.nonrep@otu_table)
SampleFilter1 <- rownames(pilot_nontree_m.nonrep@otu_table)[sample_sums(pilot_nontree_m.nonrep@otu_table)<300]
pilot_dat <- subset_samples(pilot_nontree_m.nonrep,!sample_names(pilot_nontree_m.nonrep) %in% SampleFilter1)
# save(pilot_dat, file = "pilot_dat.RDATA")

### Produce phyloseq after rarefaction
set.seed(300)
pilot_rarefied <- rarefy_even_depth(pilot_dat, sample.size=min(sample_sums(pilot_dat)), rngseed=300)
set.seed(1); .Random.seed
pilot_dat_rarefied <- prune_taxa(taxa_sums(pilot_rarefied)>0,pilot_rarefied)
# save(pilot_dat_rarefied, file = "pilot_dat_rarefied.RDATA")

##############################################################################
########################## DNA sequence analysis #############################
##############################################################################

# load("pilot_dat.RDATA") # not rarefied
load("pilot_dat_rarefied.RDATA")
pilot_dat <- pilot_dat_rarefied
pilot_dat.ReAb <- transform_sample_counts(pilot_dat, function(OTU) OTU*100/sum(OTU))
# write.csv(pilot_dat.ReAb@tax_table, "tax_table.csv")

### Correlation between read depth and weight of gut contents
pilot_study_metadata$ID <- gsub(pattern = "\\.", replacement = "-", x = pilot_study_metadata$ID)
pilot_study_metadata_merge <- merge(read_depth, pilot_study_metadata, by.x="fish_S_ID", by.y="ID")
cor.test(pilot_study_metadata_merge$read_depth, pilot_study_metadata_merge$Gut_contents_weight_g)
plot(pilot_study_metadata_merge$read_depth, pilot_study_metadata_merge$Gut_contents_weight_g)
summary(lm(read_depth ~ Gut_contents_weight_g, pilot_study_metadata_merge))
ggplot(pilot_study_metadata_merge, aes(x=Gut_contents_weight_g, y=read_depth, color=fish_species)) +
  geom_point(size=2) +
  geom_smooth(method="lm", se=FALSE, color="lightgray") +
  scale_color_manual(name="Species", values=species.cols) +
  ylab("Sample Read Depth") + xlab("Gut content weight (g)") +
  theme_bw() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.text = element_text(size=12))

### Phylum-level breakdown
# All samples
pilot_dat.ReAb.all <- prune_taxa(taxa_sums(pilot_dat.ReAb)>0, pilot_dat.ReAb)
names_unknown <- taxa_names(pilot_dat.ReAb.all@tax_table)[is.na(pilot_dat.ReAb.all@tax_table[,"Phylum"])]
pilot_dat.ReAb.all@tax_table[is.na(pilot_dat.ReAb.all@tax_table[,"Phylum"]),"Phylum"] <- c("unclassified Phyla")
pilot_dat.ReAb.all.phylum <- tax_glom(pilot_dat.ReAb.all, NArm=FALSE, taxrank="Phylum")

my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(10), rainbow(10), "azure3")
base.phylum <- pilot_dat.ReAb.all.phylum %>% psmelt()
base.phylum$Phylum <- factor(base.phylum$Phylum, levels=unique(base.phylum$Phylum))
base.phylum$fh_fz <- factor(base.phylum$fh_fz, levels=c("fresh", "frozen"), labels=c("Fresh", "Frozen"))
base.phylum$washed <- factor(base.phylum$washed, levels=c("washed", "unwashed"), labels=c("Washed", "Unwashed"))
base.phylum$fish_species <- factor(base.phylum$fish_species, levels=c("C_anomalum", "E_caeruleum"), labels=c("Herbivore", "Insectivore"))
plot2 <- ggplot(base.phylum, aes(x=fish_S_ID, y=Abundance, fill=Phylum)) + 
  geom_col(width=0.9) +
  facet_wrap(~fish_species + fh_fz + washed, scales="free_x", nrow=1, ncol=6) +
  scale_fill_manual(name="Phylum", values=my.cols) +
  scale_color_manual(name="Phylum", values=my.cols) +
  ylab("Relative Abundance (%)") + xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=8))

plot2
# ggsave("Figures/RelAbund_barplot_all_samples.pdf", plot=plot2)

### Closer look at taxonomic proportions (Phyla) per group
pie_chart_phylum <- function(dat){
  taxatable <- as.data.frame(dat@tax_table[,"Phylum"])
  taxatable$OTUname <- rownames(taxatable)
  otutable <- as.data.frame(t(as.data.frame(dat@otu_table)))
  otutable$otuname <- rownames(otutable)
  taxotutable <- merge(taxatable, otutable, by.x="OTUname", by.y="otuname")
  taxotutable <- taxotutable[,c(2:ncol(taxotutable))]
  Pgroup <- taxotutable$Phylum
  ReAb <- rowsum(taxotutable[,c(2:ncol(taxotutable))],Pgroup)
  ReAb$mReAb <- rowMeans(ReAb)
  ReAb$Phylum <- rownames(ReAb)
  return(ReAb)
}

# All samples
pie_all <- pie_chart_phylum(pilot_dat.ReAb.all.phylum)

# C_amonalum, fresh vs. frozen
dat1 <- subset_samples(pilot_dat.ReAb.all.phylum, fish_species=="C_anomalum")
dat2 <- subset_samples(dat1, fh_fz=="fresh")
pie_CA_fh <- pie_chart_phylum(dat2)
dat3 <- subset_samples(dat1, fh_fz=="frozen")
pie_CA_fz <- pie_chart_phylum(dat3)

# E_caeruleum, fresh vs. frozen
dat1 <- subset_samples(pilot_dat.ReAb.all.phylum, fish_species=="E_caeruleum")
dat2 <- subset_samples(dat1, fh_fz=="fresh")
pie_EC_fh <- pie_chart_phylum(dat2)
dat3 <- subset_samples(dat1, fh_fz=="frozen")
pie_EC_fz <- pie_chart_phylum(dat3)

### Pie charts of major Phyla
pie_plot <- merge(pie_CA_fh[,c("Phylum","mReAb")], pie_CA_fz[, c("Phylum","mReAb")], by="Phylum", all=TRUE)
names(pie_plot) <- c("Phylum","mReAb.ca.fh","mReAb.ca.fz")
pie_plot <- merge(pie_plot,pie_EC_fh[, c("Phylum","mReAb")],by="Phylum", all=TRUE)
pie_plot <- merge(pie_plot, pie_EC_fz[, c("Phylum","mReAb")], by="Phylum", all=TRUE)
names(pie_plot) <- c("Phylum","mReAb.ca.fh","mReAb.ca.fz","mReAb.ec.fh","mReAb.ec.fz")
pie_plot$MREAB <- pie_all$mReAb
pie_plot_ord <- pie_plot[order(pie_plot$MREAB, decreasing=TRUE),]
pie_plot_ord$Phylum <- factor(pie_plot_ord$Phylum,levels=c(unique(as.character(pie_plot_ord$Phylum))))
pie1 <- ggplot(pie_plot_ord,aes(x="", y=mReAb.ca.fh, fill=Phylum)) + 
  geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + 
  labs(title="C_anomalum, fresh", size=12) + 
  theme_void() + scale_fill_manual(values=c(my.cols))
pie2 <- ggplot(pie_plot_ord,aes(x="", y=mReAb.ca.fz, fill=Phylum)) + 
  geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + 
  labs(title="C_anomalum, frozen", size=12) + 
  theme_void() + scale_fill_manual(values=c(my.cols))
pie3 <- ggplot(pie_plot_ord,aes(x="", y=mReAb.ec.fh, fill=Phylum)) + 
  geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + 
  labs(title="E_caeruleum, fresh", size=12) + 
  theme_void() + scale_fill_manual(values=c(my.cols))
pie4 <- ggplot(pie_plot_ord,aes(x="", y=mReAb.ec.fz, fill=Phylum)) + 
  geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + 
  labs(title="E_caeruleum, frozen", size=12) + 
  theme_void() + scale_fill_manual(values=c(my.cols))
pie_plots_all <- ggarrange(pie1, pie2, pie3, pie4, legend='none', nrow=1, ncol=4, align="v")
pie_plots_all
# ggsave("Figures/RelAbund_piechart.pdf", plot=pie_plots_all)

### Diversity stats per sample
alpha_ASV_meta <- data.frame(Shannon_ASV=diversity(pilot_dat@otu_table,index="shannon"),
                             Observed_ASV=specnumber(pilot_dat@otu_table),
                             Pielou_ASV=diversity(pilot_dat@otu_table,index="shannon")/log(specnumber(pilot_dat@otu_table)),
                             Simpson_ASV=diversity(pilot_dat@otu_table,index="simp"),
                             Sample_ID=sample_data(pilot_dat)$fish_S_ID,
                             fh_fz=sample_data(pilot_dat)$fh_fz,
                             fish_species=sample_data(pilot_dat)$fish_species,
                             washed=sample_data(pilot_dat)$washed)
alpha_ASV_meta$fh_fz <- factor(alpha_ASV_meta$fh_fz, levels=c("fresh", "frozen"), labels=c("Fresh", "Frozen"))
alpha_ASV_meta$washed <- factor(alpha_ASV_meta$washed, levels=c("washed", "unwashed"), labels=c("Washed", "Unwashed"))
alpha_ASV_meta$fish_species <- factor(alpha_ASV_meta$fish_species, levels=c("C_anomalum", "E_caeruleum"), labels=c("Herbivore", "Insectivore"))

### Plot ASV richness: washed vs. unwashed (E. caeruleum only)
washed_df <- alpha_ASV_meta[alpha_ASV_meta$fish_species=="Insectivore",]
ASV_Richness <- ggplot(washed_df, aes(x=washed, y=Observed_ASV, fill=washed)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  facet_wrap(~fish_species, scales="free_x", nrow=1, ncol=2) +
  scale_fill_manual(values=c("#99CC99","#56B4E9")) +
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Washed w/ PBS", y="ASV richness", size=12) + 
  theme(axis.title.x=element_blank(), axis.text=element_text(size=12)) + theme_bw()
Shannon_diversity <- ggplot(washed_df, aes(x=washed, y=Shannon_ASV, fill=washed)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  facet_wrap(~fish_species, scales="free_x", nrow=1, ncol=2) +
  scale_fill_manual(values=c("#99CC99","#56B4E9")) +
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Washed w/ PBS", y="Shannon diversity", size=12) + 
  theme(axis.title.x=element_blank(), axis.text=element_text(size=12)) + theme_bw()
Pielou_Evenness <- ggplot(washed_df, aes(x=washed, y=Pielou_ASV, fill=washed)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  facet_wrap(~fish_species, scales="free_x", nrow=1, ncol=2) +
  scale_fill_manual(values=c("#99CC99","#56B4E9")) +
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Washed w/ PBS", y="Pielou's evenness", size=12) + 
  theme(axis.title.x=element_blank(), axis.text=element_text(size=12)) + theme_bw()
plot3 <- ggarrange(ASV_Richness, Shannon_diversity, Pielou_Evenness, common.legend=FALSE,
                   nrow=1, ncol=3, legend="none", align="h")
plot3
# ggsave("Figures/div_metrics_washed.pdf", plot=plot3, width=9.0, height=3.0)

### ANOVA to test differences in alpha diversity among washed and unwashed samples
richness.aov <- aov(Observed_ASV ~ washed, data=washed_df)
summary(richness.aov)
shannon.aov <- aov(Shannon_ASV ~ washed, data=washed_df)
summary(shannon.aov)
evenness.aov <- aov(Pielou_ASV ~ washed, data=washed_df)
summary(evenness.aov)
t.test(washed_df[washed_df$washed=="Washed",]$Observed_ASV, washed_df[washed_df$washed=="Unwashed",]$Observed_ASV)
t.test(washed_df[washed_df$washed=="Washed",]$Shannon_ASV, washed_df[washed_df$washed=="Unwashed",]$Shannon_ASV)
t.test(washed_df[washed_df$washed=="Washed",]$Pielou_ASV, washed_df[washed_df$washed=="Unwashed",]$Pielou_ASV)
# no significant differences 

alpha_ASV_meta <- alpha_ASV_meta[alpha_ASV_meta$washed=="washed",]
  
### Plot ASV richness: fresh vs. frozen by species
ASV_Richness <- ggplot(alpha_ASV_meta, aes(x=fh_fz, y=Observed_ASV, fill=fh_fz)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  facet_wrap(~fish_species, scales="free_x", nrow=1, ncol=2) +
  scale_fill_manual(values=c("#99CC99","#56B4E9")) +
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Fresh vs. frozen fish", y="ASV Richness", size=12) + 
  theme(axis.title.x=element_blank(), axis.text=element_text(size=12)) + theme_bw()
Shannon_diversity <- ggplot(alpha_ASV_meta, aes(x=fh_fz, y=Shannon_ASV, fill=fh_fz)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  facet_wrap(~fish_species, scales="free_x", nrow=1, ncol=2) +
  scale_fill_manual(values=c("#99CC99","#56B4E9")) +
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Fresh vs. frozen fish", y="Shannon diversity", size=12) + 
  theme(axis.title.x=element_blank(), axis.text=element_text(size=12)) + theme_bw()
Pielou_Evenness <- ggplot(alpha_ASV_meta, aes(x=fh_fz, y=Pielou_ASV, fill=fh_fz)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  facet_wrap(~fish_species, scales="free_x", nrow=1, ncol=2) +
  scale_fill_manual(values=c("#99CC99","#56B4E9")) +
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Fresh vs. frozen fish", y="Pielou's Evenness", size=12) + 
  theme(axis.title.x=element_blank(), axis.text=element_text(size=12)) + theme_bw()
plot4 <- ggarrange(ASV_Richness, Shannon_diversity, Pielou_Evenness, common.legend=FALSE,
                   nrow=1, ncol=3, legend="none", align="h")
plot4
# ggsave("Figures/div_metrics_fresh_frozen.pdf", plot=plot4, width=9.0, height=3.0)

### ANOVA to test differences in alpha diversity among fresh and frozen samples
richness.aov2 <- aov(Observed_ASV ~ fh_fz + fish_species, data=alpha_ASV_meta)
summary(richness.aov2)
shannon.aov2 <- aov(Shannon_ASV ~ fh_fz + fish_species, data=alpha_ASV_meta)
summary(shannon.aov2)
evenness.aov2 <- aov(Pielou_ASV ~ fh_fz + fish_species, data=alpha_ASV_meta)
summary(evenness.aov2)

### t-test of fresh and frozen samples when separated by species 
CA <- alpha_ASV_meta[alpha_ASV_meta$fish_species=="Herbivore",]
EC <- alpha_ASV_meta[alpha_ASV_meta$fish_species=="Insectivore",]
t.test(CA[CA$fh_fz=="Fresh",]$Observed_ASV, CA[CA$fh_fz=="Frozen",]$Observed_ASV)
t.test(CA[CA$fh_fz=="Fresh",]$Shannon_ASV, CA[CA$fh_fz=="Frozen",]$Shannon_ASV)
t.test(CA[CA$fh_fz=="Fresh",]$Pielou_ASV, CA[CA$fh_fz=="Frozen",]$Pielou_ASV)
t.test(EC[EC$fh_fz=="Fresh",]$Observed_ASV, EC[EC$fh_fz=="Frozen",]$Observed_ASV)
t.test(EC[EC$fh_fz=="Fresh",]$Shannon_ASV, EC[EC$fh_fz=="Frozen",]$Shannon_ASV)
t.test(EC[EC$fh_fz=="Fresh",]$Pielou_ASV, EC[EC$fh_fz=="Frozen",]$Pielou_ASV)
# no significant differences 

### PCOA of fresh vs. frozen samples
ord.PCoA.bray <- ordinate(pilot_dat.ReAb.all, method="PCoA", distance="bray")
pilot_dat.ReAb.all@sam_data$combined <- with(pilot_dat.ReAb.all@sam_data, paste0(fish_species, "_", fh_fz))
pcoa.cols <- c("#800080", "#800080", "#FDB81A", "#FDB81A")
plot5 <- plot_ordination(pilot_dat.ReAb.all, ord.PCoA.bray, type="samples", color="combined", shape="combined") +
  stat_ellipse(type="norm") +
  geom_point(size=6, stroke=1) + 
  scale_color_manual(values=pcoa.cols) +
  scale_shape_manual(values=c(1, 6, 1, 6)) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12)) +
  theme_bw()
plot5$layers <- plot5$layers[-1]
plot5
# ggsave("Figures/pcoa_all.pdf", plot=plot5, width=7, height=5)

plot_ordination(pilot_dat.ReAb.all, ord.PCoA.bray, type="samples", color="washed") +
  stat_ellipse(type="norm") +
  geom_point(size=6, stroke=1) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12)) +
  theme_bw()

### PERMANOVA
meta <- microbiome::meta(pilot_dat.ReAb.all)
abund <- microbiome::abundances(pilot_dat.ReAb.all)
meta_CA <- meta[meta$fish_species=="C_anomalum",]
abund_CA <- abund[,colnames(abund) %in% meta_CA$fish_S_ID]
meta_EC <- meta[meta$fish_species=="E_caeruleum",]
abund_EC <- abund[,colnames(abund) %in% meta_EC$fish_S_ID]
adonis2(t(abund_CA) ~ fh_fz, data=meta_CA, permutations=999, method="bray")
adonis2(t(abund_EC) ~ fh_fz, data=meta_EC, permutations=999, method="bray")

### PERMDISP
dist_adund <- phyloseq::distance(pilot_dat.ReAb.all, method="bray")
dist_adund_CA <- phyloseq::distance(subset_samples(pilot_dat.ReAb.all,fish_species=="C_anomalum"), method="bray")
dist_adund_EC <- phyloseq::distance(subset_samples(pilot_dat.ReAb.all,fish_species=="E_caeruleum"), method="bray")
beta_CA <- betadisper(dist_adund_CA, meta_CA$fh_fz)
permutest(beta_CA)
beta_EC <- betadisper(dist_adund_EC, meta_EC$fh_fz)
permutest(beta_EC)

### CCA 
pilot.cca <- ordinate(pilot_dat, "CCA")
p <- plot_ordination(pilot_dat, pilot.cca,
                     type = "split", shape = "fish_species", 
                     color = "Phylum", label = "fh_fz")
p + geom_point(size = 4)

### Fresh vs. frozen LDA, Family level
fh.fz.cols <- c("salmon", "lightblue")
pilot_dat_t <- phyloseq(otu_table(t(otu_table(pilot_dat))), tax_table(pilot_dat), sample_data(pilot_dat))
LL <- 1:nrow(pilot_dat_t@tax_table)
for(i in LL){
  if(is.na(pilot_dat_t@tax_table[i,2])){
    pilot_dat_t@tax_table[i,5] <- paste0("unclassified ","Bacteria"," lineage")
    pilot_dat_t@tax_table[i,4] <- c("unclassified Pylum")
    pilot_dat_t@tax_table[i,3] <- c("unclassified Phylum")
  }else{
    if(is.na(pilot_dat_t@tax_table[i,3])){
      pilot_dat_t@tax_table[i,5] <-  paste0("unclassified ",pilot_dat_t@tax_table[i,2]," lineage")
      pilot_dat_t@tax_table[i,4] <-  c("unclassified Class")
      
    }else{
      if(is.na(pilot_dat_t@tax_table[i,4])){
        pilot_dat_t@tax_table[i,5] <-  paste0("unclassified ",pilot_dat_t@tax_table[i,3]," lineage")
        
      }else{
        if(is.na(pilot_dat_t@tax_table[i,5])){
          pilot_dat_t@tax_table[i,5] <-  paste0("unclassified ",pilot_dat_t@tax_table[i,4]," lineage")
        }
      }
    }
  }
}
pilot_dat_t.F <- tax_glom(pilot_dat_t, NArm=FALSE, taxrank = "Family")
fh_fz_LDA <- phyloseq_to_deseq2(pilot_dat_t.F, ~fh_fz)
set.seed(9876)
res <- lefser(fh_fz_LDA, groupCol="fh_fz", lda.threshold=0.5)
res$Names <- pilot_dat_t.F@tax_table[res$Names,5]
head(res)
paste(res[res$scores>0,]$Names,collapse = ", ")
paste(res[res$scores<0,]$Names,collapse = ", ")
groupf <- colData(fh_fz_LDA)[["fh_fz"]]
groupf <- as.factor(groupf)
groupsf <- levels(groupf)
.numeric01 <- function(x) {
  x <- as.factor(x)
  uvals <- levels(x)
  ifelse(x == uvals[1L], 0L, 1L)
}
group <- .numeric01(groupf)
lefse_plot1 <- lefserPlot(res) + 
  labs(fill="Fresh vs. frozen", title="Family") + 
  scale_fill_manual(values = fh.fz.cols, labels=c("Fresh","Frozen"))
lefse_plot1

### Fresh vs. frozen LDA, Genus level
pilot_dat_t.G <- tax_glom(pilot_dat_t, NArm=TRUE, taxrank="Genus")
fh_fz_LDA <- phyloseq_to_deseq2(pilot_dat_t.G , ~fh_fz)
res <- lefser(fh_fz_LDA, groupCol="fh_fz", lda.threshold=0.5)
res$Names <- pilot_dat_t.G@tax_table[res$Names,6]
head(res)
paste(res[res$scores>0,]$Names,collapse = ", ")
paste(res[res$scores<0,]$Names,collapse = ", ")
groupf <- colData(fh_fz_LDA)[["fh_fz"]]
groupf <- as.factor(groupf)
groupsf <- levels(groupf)
group <- .numeric01(groupf)
lefse_plot2 <- lefserPlot(res) + 
  labs(fill="Fresh vs. frozen", title="Genus") + 
  scale_fill_manual(values=fh.fz.cols, labels=c("Fresh","Frozen"))
lefse_plot2

### Herbivore vs. insectivore LDA, Family level
spec.cols <- c("#800080", "#FDB81A")
spec_LDA <- phyloseq_to_deseq2(pilot_dat_t.F, ~fish_species)
res <- lefser(spec_LDA, groupCol="fish_species", lda.threshold=1.6)
res$Names <- pilot_dat_t.F@tax_table[res$Names,5]
head(res)
paste(res[res$scores>0,]$Names,collapse = ", ")
paste(res[res$scores<0,]$Names,collapse = ", ")
groupf <- colData(spec_LDA)[["fish_species"]]
groupf <- as.factor(groupf)
groupsf <- levels(groupf)
.numeric01 <- function(x) {
  x <- as.factor(x)
  uvals <- levels(x)
  ifelse(x == uvals[1L], 0L, 1L)
}
group <- .numeric01(groupf)
lefse_plot3 <- lefserPlot(res) + 
  labs(title="Family") + 
  scale_fill_manual(values = spec.cols, labels=c("Herbivore","Insectivore")) +
  theme_bw()
lefse_plot3
# ggsave("Figures/lefse_species.pdf", plot=lefse_plot3)

### Herbivore vs. insectivore LDA, Genus level
pilot_dat_t.G <- tax_glom(pilot_dat_t, NArm=TRUE, taxrank="Genus")
spec_LDA <- phyloseq_to_deseq2(pilot_dat_t.G , ~fish_species)
res <- lefser(spec_LDA, groupCol="fish_species", lda.threshold=1.5)
res$Names <- pilot_dat_t.G@tax_table[res$Names,6]
head(res)
groupf <- colData(spec_LDA)[["fish_species"]]
groupf <- as.factor(groupf)
groupsf <- levels(groupf)
group <- .numeric01(groupf)
lefse_plot4 <- lefserPlot(res) + 
  scale_fill_manual(values=spec.cols, labels=c("Herbivore","Insectivore"))
lefse_plot4

