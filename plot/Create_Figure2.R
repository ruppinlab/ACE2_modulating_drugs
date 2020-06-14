# This script is used to generate Figure 2-main text and Figure EV3
require(ggplot2)
require(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

################################################
# Panel 2A
################################################
new_collection=readRDS('../data/Figure2A_collected.results.with.old.ones.RDS') #! please fix: file not there
new_collection_df=data.frame(padj=new_collection$padj, 
                             logfc=new_collection$log.fc,
                             drugName=new_collection$drug)
# Setting a max cap on the LFC change
new_collection_df$logfc[new_collection_df$logfc>999]=5
new_collection_df$point_color=''
new_collection_df[new_collection_df$logfc > 0 & -log(new_collection_df$padj, 10)>1.3,]$point_color='positive'
new_collection_df[new_collection_df$logfc < -0 & -log(new_collection_df$padj, 10)>1.3,]$point_color='negative'
Figure2A<-ggplot(new_collection_df, aes(y= -log(padj, 10), x=logfc, color=point_color))+
  geom_point(size=4)+
  theme_bw(base_size = 20)+
  labs(y='-log10(P)', x='ACE2 expression logFC')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1.5, color='blue')+
  geom_label_repel(data=subset(new_collection_df, -log(padj, 10)>1.3),
                   aes(label=paste(drugName)), size=6)+
  ggtitle('Lung in vivo & in vitro from GEO')+
  theme(legend.position = 'none')+
  annotate("text", x = -1.5, y = 1.4, label = "P < 0.05", color='blue', size=5)+
  scale_color_manual(values = c('black', '#F8766D', '#00BFC4'))
################################################
# Panel 2B
################################################
df_completev2=readRDS('../data/figure2B_Data.RDS')
tiff('Figure2C_nonCancer.tiff',width = 500, height=700)
Figure2B<-ggplot(df_completev2[df_completev2$whether_cancer=='non-cancer',],
                 aes(y= ACE2, fill=Group,x=drug, label=pval))+
  geom_boxplot()+
  geom_text()+
  theme_bw(base_size = 20)+
  theme(legend.position = 'top',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y='ACE2 Expression')+
  facet_wrap(~drug_andTissue+GEO, shrink = T, scales = 'free', nrow = 2)
dev.off()
################################################
# Panel 2C
################################################
kidney_collection=readRDS('../data/Figure2C_kidney.collected.results.RDS') #! please fix: file not there
kidney_collection_df=data.frame(padj=kidney_collection$padj, 
                                logfc=kidney_collection$log.fc,
                                drugName=kidney_collection$drug)
kidney_collection_df$drugName=as.character(kidney_collection_df$drugName)
kidney_collection_df[grep('cisplatin',kidney_collection_df$drug),]
kidney_collection_df$drugName[grep('cisplatin',kidney_collection_df$drugName)][1:2]=
  c('cisplatin (Renal\n cortex, rat)', 'cisplatin (Whole\n kidney, mice)')
kidney_collection_df$drugName=factor(kidney_collection_df$drugName)
kidney_collection_df$logfc[kidney_collection_df$logfc>9]=5
kidney_collection_df$point_color=''
kidney_collection_df[kidney_collection_df$logfc > 0 & -log(kidney_collection_df$padj, 10)>1.3,]$point_color='positive'
kidney_collection_df[kidney_collection_df$logfc < -0 & -log(kidney_collection_df$padj, 10)>1.3,]$point_color='negative'
Figure2C<-ggplot(kidney_collection_df, aes(y= -log(padj, 10), x=logfc, color=point_color))+
  geom_point(size=4)+
  theme_bw(base_size = 20)+
  labs(y='-log10(P)', x='ACE2 expression logFC')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1.5, color='blue')+
  geom_label_repel(data=subset(kidney_collection_df, -log(padj, 10)>1.3),
                   aes(label=paste(drugName)), size=6)+
  ggtitle('Kidney in vivo & in vitro from GEO')+
  theme(legend.position = 'none')+
  annotate("text", x = -0.50, y = 1.4, label = "P < 0.05", color='blue', size=5)+
  scale_color_manual(values = c('black', '#F8766D', '#00BFC4'))
################################################
# Plot
################################################
df_complete_kidney=readRDS('../data/figure2B_Data.RDS')
tiff('Figure2D_Kidney.tiff',width = 600, height=900)
Figure2D<-ggplot(df_complete_kidney, aes(y= ACE2, fill=Group,x=drug, label=pval))+
  geom_boxplot()+
  geom_text()+
  theme_bw(base_size = 20)+
  theme(legend.position = 'top',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y='ACE2 Expression')+
  facet_wrap(~drug_andTissue+GEO, shrink = T, scales = 'free', nrow = 3)
dev.off()
#################
# Figure EV3
#################
tiff('EV3_Cancer.tiff',width = 800, height=800)
ggplot(df_completev2[df_completev2$whether_cancer=='cancer',],
       aes(y= ACE2, fill=Group,x=drug, label=pval))+
  geom_boxplot()+
  geom_text()+
  theme_bw(base_size = 20)+
  theme(legend.position = 'top',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y='ACE2 Expression')+
  facet_wrap(~drug_andTissue+GEO, shrink = T, 
             scales = 'free', nrow = 3)
dev.off()
