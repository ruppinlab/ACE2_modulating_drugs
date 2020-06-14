# Figure 1 without Panel D
load('../data/Figure1_new.results.RData')
require(ggplot2)
require(ggrepel)
require(readxl)
require(fgsea)
###########################################################################
# Fig-A
###########################################################################
de.res.hyp.carc$point_color=''
de.res.hyp.carc[de.res.hyp.carc$log.fc > 0 & -log(de.res.hyp.carc$pval, 10)>1,]$point_color='positive'
de.res.hyp.carc[de.res.hyp.carc$log.fc < -0 & -log(de.res.hyp.carc$pval, 10)>1,]$point_color='negative'
figA <- ggplot(de.res.hyp.carc, aes(y= -log(pval, 10), x=log.fc, color=point_color))+
  geom_point(size=4)+
  theme_bw(base_size = 20)+
  labs(y='-log10(P)', x='ACE2 expression logFC')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1, color='blue')+
  geom_label_repel(data=subset(de.res.hyp.carc, -log(pval, 10)>1.5),
                   aes(label=paste(pert_iname)), size=6)+
  ggtitle('Approved anti-hypertension drugs')+
  theme(legend.position = 'none')+
  annotate("text", x = -0.5, y = 1.35, label = "P<0.05", color='blue', size=5)+
  scale_color_manual(values = c('black', '#F8766D'))
###########################################################################
# Fig-B
###########################################################################
gsea.res.hyp.carc=get(load('../data/Figure_1B_gsea.res1.RData'))
gsea.res.hyp.carc$Enriched_in='negative modifier'
gsea.res.hyp.carc$Enriched_in[gsea.res.hyp.carc$ES>0]='positive modifier'
gsea.res.hyp.carc=gsea.res.hyp.carc[order(gsea.res.hyp.carc$pval),]
gsea.res.hyp.carc$Enriched_in=factor(gsea.res.hyp.carc$Enriched_in, labels=c('ACE2 down-regulator',
                                                                             'ACE2 up-regulator'))
fig1B<-ggplot(gsea.res.hyp.carc, aes(x=reorder(pathway, -pval, FUN = median),
                                     y= -log(pval, 10), fill=Enriched_in))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  theme_bw(base_size = 20)+
  theme(axis.text.y = element_text(angle = 0), legend.position = 'top', legend.title = element_blank(),
        legend.justification='left', legend.direction='vertical')+
  labs(x='', y='-log10(P)')+
  ggtitle('Enrichment by MOA')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1.5, color='blue')
###########################################################################
# Fig-C
###########################################################################
de.res.app.carc$point_color=''
de.res.app.carc[de.res.app.carc$log.fc > 0 & -log(de.res.app.carc$padj, 10)>1,]$point_color='positive'
de.res.app.carc[de.res.app.carc$log.fc < -0 & -log(de.res.app.carc$padj, 10)>1,]$point_color='negative'
figC <- ggplot(de.res.app.carc, aes(y= -log(pval, 10), x=log.fc, color=point_color))+
  geom_point(size=4)+
  theme_bw(base_size = 20)+
  labs(y='-log10(P)', x='ACE2 expression logFC')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1.5, color='blue')+
  geom_label_repel(data=subset(de.res.app.carc, -log(padj, 10)>1),
                   aes(label=paste(pert_iname)), size=6)+
  ggtitle('All approved Drugs')+
  theme(legend.position = 'none')+
  annotate("text", x = -0.50, y = 1.4, label = "P < 0.05", color='blue', size=5)+
  scale_color_manual(values = c('black', '#F8766D', '#00BFC4'))
###########################################################################
# Fig-E
###########################################################################
gsea.res_all_subset=read_xlsx('Figure_1E_topPathways_allClinicalApproveddrugsACE2.xlsx') #! please fix: file not there; maybe replace with a RData or RDS file
colnames(gsea.res_all_subset)=colnames(gsea.res_all)
gsea.res_all_subset$Enriched_in='negative modifier'
gsea.res_all_subset$Enriched_in[gsea.res_all_subset$ES>0]='positive modifier'
figE<-ggplot(gsea.res_all_subset[1:10,], aes(x=reorder(pathway, -pval), y= -log(pval, 10), fill=Enriched_in))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  theme_bw(base_size = 20)+
  theme(axis.text.y = element_text(angle = 0), legend.position = 'none', legend.title = element_blank())+
  labs(x='', y='-log10(P)')+
  ggtitle('Enrichment by MOA')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1.5, color='blue')
###########################################################################
# Fig-F
###########################################################################
load('../data/gsea.res_Figure1F.RData')
gsea.res=gsea.res[nchar(pathway)==3 & size>10][, padj:=p.adjust(pval, "BH")][pval<0.05][order(pval)]
gsea.res=gsea.res[order(gsea.res$pval),]
gsea.res$Enriched_in='negative modifier'
gsea.res$Enriched_in[gsea.res$ES>0]='positive modifier'
gsea.res$name=gsub('AND','&',gsea.res$name)
figF<-ggplot(gsea.res, aes(x=reorder(tolower(name), -pval), y= -log(pval, 10), fill=Enriched_in))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  theme_bw(base_size = 20)+
  theme(axis.text.y = element_text(angle = 0),
        legend.position = 'none', legend.title = element_blank())+
  labs(x='Indication', y='-log10(P)')+
  
  ggtitle('Enrichment by Indication')+
  geom_hline(yintercept = 1.3, linetype='dashed', size=1.5, color='blue')
###########################################################################
# "One func to combine them all" - LOTR
###########################################################################
pdf('FIG1_Compv10_wdLabels.pdf', width = 18, height=24)
plot_grid(figA,  fig1B, figC,figE, figF, nrow=3,  labels = c('A', 'B', 'C', 'D', 'E'), label_size = 32)
dev.off()
