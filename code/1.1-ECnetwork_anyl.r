library(openxlsx)
library(dplyr)
#library(tidyverse)
library(igraph)
library(reshape2)
library(gplots)
library(cowplot)
library(ggrepel)
library(ggExtra)
library(dplyr)
library(forcats)

setwd('G:/My Drive/00_professional/projects/Coronavirus/Microbiome ~ Viral Load Study (Cameron)/')

reload=F
tpm_norm = F
makeplt=T
genescore_normalization = F ### important: changes pathway activation results (if T, normalize by genus, if F, skewed by total counts (appropriate))
taxas = c('genus','sp')[2] ### select taxa

if(reload){
  # load data
  if(!file.exists('microbiome-mapping/feature-metadata-and-ranks.csv')){
    stop('Microbiome samples contain protected patient data. Please request these files through the FINRISK secure data portal')
  }
  dat = read.csv('microbiome-mapping/feature-metadata-and-ranks.csv')
  dat$Cross.reference..CAZy. = gsub(';$','',dat$Cross.reference..CAZy.)
  dat = separate_rows(dat,Cross.reference..CAZy.,sep=';')
  
  # load id mapping
  ec = read.csv('EC2CAZy.csv')
  
  # integrate mapping with data
  tmp_dat = merge(ec,dat,by.x='CAZy.ID',by.y='Cross.reference..CAZy.',all.y=T)
  
  write.csv(tmp_dat,file='microbiome-mapping/feature-metadata-and-ranks.complete_EC_CAZyID.csv',row.names = F)
  
  
  # load ec network
  ecnet = read.xlsx('code/glycogene-clusters/CAZy/cazyEC.xlsx',sheet = 3)
  ecnet_tmp = unique(melt(ecnet,id.vars = c('Target'))[,c(1,3)])
  head(ecnet)
  
  save(ecnet,ecnet_tmp,dat,tmp_dat,file='pathway-activation/ecnet.rda')
}else{
  load('pathway-activation/ecnet.rda')
}

# count normalization

rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}




for(taxa in taxas){
  
  # integrate data into the ec network
  if(tpm_norm){
    tmp_dat$log_counts = log(tpm(tmp_dat$total_counts,tmp_dat$Length))
  }else{
    tmp_dat$log_counts = log(tmp_dat$total_counts+1)  ## changed: log(counts) before all calculations
  }
  tmp_dat$taxa = tmp_dat[[taxa]]
  reactions = na.omit(merge(ecnet_tmp,unique(tmp_dat[,c('total_counts','EC','taxa','log_counts','genus','sp')]),by.x='value',by.y='EC',all.x=T))
  reactions_tmp = reactions
  
  
  ########
  pos <- position_jitter(width = 0.2,height=.01, seed = 1)
  
  top = unique(reactions$taxa[order(-reactions$log_counts)][1:150])
  ggplot(droplevels(unique(reactions[,c('Target','log_counts','sp')])),aes(x=Target,y=log_counts))+geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color=ifelse(sp%in%top,as.character(sp),NA),
                   alpha=ifelse(sp%in%top,1,.1)),
               position=pos)
  
  if(makeplt){
    ########
    sz=15
    ### variance stabilization
    g1=ggplot(reactions[reactions$genus %in% reactions$genus[order(-reactions$total_counts)[1:500]],],
           aes(x=genus,y=total_counts))+geom_boxplot()+ylab('counts by EC')+xlab('genus (top 500 species)')+
      theme_bw(base_size = sz)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
    g2=ggplot(reactions[reactions$genus %in% reactions$genus[order(-reactions$total_counts)[1:500]],],
           aes(x=genus,y=total_counts))+geom_boxplot()+scale_y_log10()+ylab('log counts by EC')+xlab('genus (top 500 species)')+
      theme_bw(base_size = sz)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
    
    tmp = reactions %>% group_by(value,taxa) %>% summarise(EC_abundance_sum = sum(exp(log_counts),na.rm=T),EC_abundance_mean = mean(exp(log_counts),na.rm=T),EC_abundance_var = var(total_counts,na.rm=T))
    g3=ggplot(tmp[tmp$EC_abundance_var>0,],aes(x=EC_abundance_mean,y=EC_abundance_var))+geom_point()+
      scale_x_log10()+scale_y_log10()+stat_smooth()+xlab('log(mean(counts)) by EC')+ylab('log(var(counts)) by EC')+theme_bw(base_size = sz)
    tmp = reactions %>% group_by(value,taxa) %>% summarise(EC_abundance_sum = sum(log_counts,na.rm=T),EC_abundance_mean = mean(log_counts,na.rm=T),EC_abundance_var = var(log_counts,na.rm=T))
    g4=ggplot(tmp[tmp$EC_abundance_var>0,],aes(x=EC_abundance_mean,y=EC_abundance_var))+geom_point()+
      scale_x_log10()+scale_y_log10()+stat_smooth()+xlab('log(mean(log counts)) by EC')+ylab('log(var(log counts)) by EC')+theme_bw(base_size = sz)
    
    g = plot_grid(g1,g2,g3,g4,labels=LETTERS[1:4])
    ggsave(g,filename='pathway-activation/QC/01-variance_stabilization.pdf',height=10,width=8)
  }  
  
  ### sum, mean and mean-threshold on log counts
  reactions = reactions %>% group_by(Target,value,taxa) %>% mutate(EC_abundance_sum = sum(log_counts,na.rm=T))
  reactions = reactions %>% group_by(Target,value,taxa) %>% mutate(EC_abundance_mean = mean(log_counts,na.rm=T))
  if(nrow(reactions)!=nrow(reactions_tmp) | nrow(reactions)<10){stop('matrix mutation has made a mistake')}
  
  if(makeplt){
      #### check distributions
    g=ggplot( reactions ,aes(y=EC_abundance_sum,x=EC_abundance_mean) )+geom_point()+
      geom_vline(aes(xintercept = median(EC_abundance_mean)),size=2,col='black')+geom_hline(aes(yintercept = median(EC_abundance_sum)),size=2,col='black')+
      geom_vline(aes(xintercept = mean(EC_abundance_mean)),size=2,col='lightblue')+geom_hline(aes(yintercept = mean(EC_abundance_sum)),size=2,col='lightblue')+
      scale_y_log10()+scale_x_log10()+theme_bw(base_size=sz)
    g=ggExtra::ggMarginal(g, type = "histogram")
    g
    ggsave(g,filename='pathway-activation/QC/02-abundance_agg_distributions.pdf',height=6,width=6)
    
    g1=ggplot(reactions,aes(x=Target,y=EC_abundance_sum))+geom_boxplot()+
      ylab('sum(log counts) by Taxa, Pathway & EC')+theme_bw(base_size = sz)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
    g2=ggplot(reactions,aes(x=Target,y=EC_abundance_mean))+geom_boxplot()+
      ylab('mean(log counts) by Taxa, Pathway & EC')+theme_bw(base_size = sz)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
    g = plot_grid(g1,g2,labels=LETTERS[1:2])
    g
    ggsave(g,filename='pathway-activation/QC/03-abundance_agg_distributions.by_pathway.pdf',height=6,width=6)
  }  
  
  reactions = reactions %>% group_by(Target,taxa) %>% 
          mutate(pathway_completeness_abundance = as.numeric(sum(EC_abundance_sum>quantile(reactions$EC_abundance_sum,.05))==length(EC_abundance_sum)))
  

  unique(reactions[reactions$Target=='Heparan Sulfate',c('genus','sp','pathway_completeness_abundance')]) %>%
    mutate(genus = fct_reorder(genus, pathway_completeness_abundance,.fun='sum')) %>%
      ggplot( aes(x=genus,fill=pathway_completeness_abundance==1))+geom_bar()
  unique(reactions[reactions$Target=='Heparan Sulfate'&reactions$genus=='Bacteroides',c('genus','sp','pathway_completeness_abundance')])
  
  #reactions = unique(reactions)
  #g1 = ggplot(reactions,aes(EC_abundance_sum,EC_abundance_mean,color=Target))+geom_point()+scale_x_log10()+theme_bw(base_size = 15)
  
  # EC_abundance_sum: sum of log_counts within each taxa for each reaction within each pathway
  # EC_abundance_mean: mean of log_counts within each taxa for each reaction within each pathway
  # pathway_completeness_abundance: 1 or 0 corresponding indicating if all reaction in a pathway are active (>global mean)
  
  
  ### calculate gene_scores: gene abundance relative to typical abundance for that EC : c*log(abundance/bounded-average)
  range = c(.25,.75)
  reactions = reactions %>% group_by(genus) %>% mutate(threshold_p = ecdf(log_counts)(mean(log_counts)))
  reactions = reactions %>% group_by(genus) %>% mutate(threshold_p_adj = ifelse(threshold_p<.5, max(range[1],threshold_p),min(range[2],threshold_p) ) )
  reactions = reactions %>% group_by(genus) %>% mutate(threshold = quantile(log_counts,threshold_p_adj) )
  if(genescore_normalization){
    reactions$gene_score = 5*log(1+(reactions$log_counts/reactions$threshold))   #
  }else{
    reactions$gene_score = reactions$log_counts
  }

  if(nrow(reactions)!=nrow(unique(reactions)) | nrow(reactions)<10){stop('matrix mutation has made a mistake')}
  
  if(makeplt){
    ##### TODO MAKE A FIGURE
    g2a <- ggplot(reactions,aes(x=EC_abundance_mean,fill=Target))+geom_density(alpha=.2)+theme_bw(base_size = 15)
    g2b <- ggplot(reactions,aes(x=gene_score,fill=Target))+geom_density(alpha=.2)+theme_bw(base_size = 15)
    
    g = plot_grid(g2a,g2b,labels=LETTERS[1:2])
    g
    ggsave(g,filename='pathway-activation/QC/04-gene_activation_distributions.by_pathway.pdf',height=4,width=6)
  }
  
  # threshold: 25th/75th percentile bounded mean within a CAZy ID; abundance[percentile(mean s.t. 25%<=mean<=75%)] (not used)
  # gene_score: constant-scaled log ratio of abundance to threshold
  
  ### calculate pathway capacity: for reactions in pathway, min reaction_score, max gene_score
  reactions = reactions %>% group_by(Target,value,taxa) %>% mutate(EC_score_max = max(gene_score,na.rm=T))
  reactions = reactions %>% group_by(Target,value,taxa) %>% mutate(EC_score_mean = mean(gene_score,na.rm=T))
  reactions = reactions %>% group_by(Target,taxa) %>% mutate(pathway_capacity = min(EC_score_max,na.rm=T))
  
  if(nrow(reactions)!=nrow(unique(reactions)) | nrow(reactions)<10){stop('matrix mutation has made a mistake')}
  
  if(makeplt){
    ### data exploration
    
    g3 = ggplot(reactions,aes(EC_abundance_mean,pathway_capacity,color=Target))+geom_point()+theme_bw(base_size = 15)
    g4 = ggplot(reactions,aes(x=Target,fill=factor(pathway_completeness_abundance),pathway_capacity))+geom_boxplot()+theme_bw(base_size = 15)+
      scale_fill_brewer(name='pathway\ncompleteness')
    
    # EC_score_max: maximum activity score for a reaction in a pathway within a taxa
    # EC_score_mean: average activity score for a reaction in a pathway within a taxa
    # pathway_capacity: estimated capacity for each pathway; min_reaction(max_gene(gene activation score))
    
    g=plot_grid(g3,g4)
    g
    ggsave(g,filename = paste0('pathway-activation/pathway_completeness-capacity.',taxa,'.pdf'),height = 20,width=20)
  
    
    rank = with(reactions, data.frame(genus=genus,sp=sp,value=pathway_capacity))
    top = unique(rank$sp[order(rank$value,decreasing = TRUE)[1:200]])
    top5k = unique(rank$sp[order(rank$value,decreasing = TRUE)[1:1000]])
    
    pos <- position_jitter(width = 0.2,height=.01, seed = 1)
    g5 = ggplot(na.omit(droplevels(unique(reactions[,c('pathway_capacity','Target','sp','genus','EC_abundance_sum')]))),
                aes(y=pathway_capacity,x=Target,label=ifelse(sp%in%top,as.character(sp),NA)))+
      geom_boxplot(outlier.shape = NA)+
      #geom_text_repel(position = pos)+
      geom_point(aes(color=ifelse(sp%in%top,as.character(sp),NA),alpha=ifelse(sp%in%top,1,.1)),position=pos)+
      guides(alpha = F)+
      labs(color = 'Species')+
      theme_bw(base_size = 15)#+theme(legend.position = "top")#+scale_color_viridis_d()
    g5
    ggsave(g5,filename = paste0('pathway-activation/pathway_completeness-capacity.',taxa,'.named.pdf'),height=10,width=15)
    
    rank = with(reactions, data.frame(genus=genus,sp=sp,pathway_capacity=pathway_capacity))
    rank$rank = order(rank$pathway_capacity)
  
    g6=ggplot(droplevels(unique(reactions[reactions$sp%in%as.character(top),c('pathway_capacity','genus','sp','Target')])),
           aes(x=Target,fill=genus))+
      geom_bar()+
      #geom_bar(stat='identity',position='dodge')+
      scale_fill_viridis_d()+theme_bw(base_size=15)
    g6
    ggsave(g6,filename = paste0('pathway-activation/pathway_completeness-capacity.',taxa,'.barplots.pdf'))
  }
  
  reactions = unique(reactions)
  
  write.csv(reactions,file=paste0('pathway-activation/pathway_completeness-capacity.',taxa,'.csv'))
  
  #### additional data exploration
  d1 = dcast(data=unique(reactions[,c('Target','value','taxa','EC_abundance_sum')]),formula= Target+value~taxa,value.var = 'EC_abundance_sum',fill=0)
  d2 = dcast(data=unique(reactions[,c('Target','taxa','pathway_completeness_abundance')]),formula= Target~taxa,value.var = 'pathway_completeness_abundance',fill=0,fun.aggregate = prod)
  
  d1b = dcast(data=unique(reactions[,c('Target','value','taxa','EC_score_max')]),formula= Target+value~taxa,value.var = 'EC_score_max',fill=0)
  d2b = dcast(data=unique(reactions[,c('Target','taxa','pathway_capacity')]),formula= Target~taxa,value.var = 'pathway_capacity',fill=0,fun.aggregate = prod)
  
  if(makeplt){
    pdf(paste0('pathway-activation/pathway_completeness.',taxa,'.heatmaps.pdf'),height=20)
    heatmap.2(log(t(data.matrix(d1[,-c(1:3)]))+1),Colv = F,trace='none',ColSideColors = rainbow(7)[factor(d1$Target)],col=hcl.colors,
              labCol = paste(d1$value,d1$Target),cexCol=1,margins=c(12,7))
    heatmap.2(t(data.matrix(d2[,-1])),trace='none',labCol = d2$Target,margins = c(10,7),cexCol = 1,col=rev(c('black','lightgrey')))
    dev.off()
    
    pdf(paste0('pathway-activation/pathway_capacity.',taxa,'.heatmaps.pdf'),height=20)
    heatmap.2(t(data.matrix(d1b[,-c(1:3)])),Colv = F,trace='none',ColSideColors = rainbow(7)[factor(d1$Target)],col=hcl.colors,
              labCol = paste(d1$value,d1$Target),cexCol=1,margins=c(12,7))
    heatmap.2(t(data.matrix(d2b[,-1])),trace='none',labCol = d2$Target,margins = c(10,7),cexCol = 1,col=hcl.colors)
    dev.off()
    
    pdf(paste0('pathway-activation/pathway_capacity.',taxa,'.heatmaps.small.pdf'),height=20)
    rank = with(reactions, data.frame(genus=genus,sp=sp,value=pathway_capacity))
    top = unique(rank$sp[order(rank$value,decreasing = TRUE)[1:1000]])
    
    d2b_reduced = t(data.matrix(d2b[,-1]))
    d2b_reduced = d2b_reduced[rownames(d2b_reduced) %in% top,]

    rownames(d2b_reduced) = factor(gsub(' .*','',rownames(d2b_reduced)))
    genus = factor(rownames(d2b_reduced))
    
    heatmap.2(t(d2b_reduced),trace='none',labRow = d2$Target,
              margins = c(10,7),cexCol = 1,col=hcl.colors,
              ColSideColors = rainbow(length(levels(genus)))[genus]
    )
    dev.off()
    
    bacs = c('Bacteroides ovatus','Bacteroides thetaiotaomicron','Flavobacterium heparinum')
    distr_cap = ecdf(reactions$pathway_capacity)
    hist(perc<-(1-distr_cap(reactions$pathway_capacity)))
    hist(p.adjust(perc,'fdr'))
    unique(cbind(as.character(reactions$sp),perc)[reactions$sp%in%bacs,])
  }
}


