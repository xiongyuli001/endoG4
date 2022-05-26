# 1.structure stablity ----------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
pqs <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed') %>% data.frame()
G4 <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed') %>% data.frame()
colnames(pqs)[5] <- 'score'
pqs$type <- 'PQS'

colnames(G4)[5] <- 'score'
G4$type <- 'G4'
df <- rbind(pqs,G4)
#计算中位数
medians <- aggregate(score ~  type, df, median)
a <- table(df$type) %>% data.frame()
colnames(a)[1] <- 'type'
colnames(a)[2] <- 'number'
p <- ggplot(data = df,aes(x = type,y = score,fill = type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white",alpha = 0.8) +
  #geom_text(data = medians, aes(label = round(score,3), y = 85),size = 8,col = 'red')+
  #geom_text(data = a,aes(label = as.character(number),y=95),size = 8,col = 'black')+
  scale_x_discrete(limits=c("PQS","G4")) +
  stat_compare_means(comparisons = list(c("PQS","G4")),label.y = 100,
                     aes(label = paste0("p = ", ..p.format..)),
                     tip.length = 0,size=5)+
  scale_fill_manual(values = c('#b2182b','#4d4d4d') ) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 18),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("Structural stability")+
  coord_cartesian(ylim=c(50,120))+
  geom_hline(yintercept = medians[2,2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/structure_stability.PQS_G4.pdf",
                  width = 4,height = 3,devices = NULL,res = 300,show = F)   
# 2.ChIPseeker -------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
library(ggplot2)
library(data.table)
pqs <- '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed'
G4 <- '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed'
peakAnno_pqs <- annotatePeak(pqs, tssRegion=c(-1500, 1500),
                             TxDb=txdb, annoDb="org.Hs.eg.db")  
peakAnno_G4 <- annotatePeak(G4, tssRegion=c(-1500, 1500),
                            TxDb=txdb, annoDb="org.Hs.eg.db")
a_pqs <- peakAnno_pqs@annoStat
a_peak <- peakAnno_G4@annoStat
a_pqs$type <- 'Non-eG4'
a_peak$type <- 'G4'
df <- rbind(a_peak,a_pqs)
rownames(df) <- 1:18
a <- vector()
a[1] <- df[4,2]+df[5,2]
a[2] <- df[6,2]+df[7,2]
a[3] <- df[13,2]+df[14,2]
a[4] <- df[15,2]+df[16,2]

df2 <- data.frame(Feature = c('Exon','Intron','Exon','Intron'),
                  Frequency =a,
                  type = c('eG4','eG4','Non-eG4','Non-eG4'))
df3 <- df[-c(4:7,13:16),]
df3$type <- if_else(df3$type == 'G4','eG4',df3$type)
df4 <- rbind(df2,df3)

p <- ggplot(df4,aes(x = Feature,y= Frequency,fill = type)) +
  geom_bar(stat = 'identity',position=position_dodge(1),width=1,color = 'white')+
  scale_fill_manual(values = c('#b2182b','#4d4d4d') ) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.x = element_text(size = 18),
        axis.title.y =element_blank(),legend.position="top")+
  ylab("Frequency (%)")+
  geom_text(aes(label=round(Frequency,1),y=Frequency+1.8), color="black", size=3.5,
            position = position_dodge(.9))+
  coord_flip(ylim = c(0,50))

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/feature.frequency.pdf",
                  width = 7.5,height = 5,devices = NULL,res = 300,show = F)  

# 3.number_of_G4------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
hg_pqs <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.result')

a <- table(hg_pqs$sum) %>% data.frame()
a <- a[-1,]
colnames(a) <- c('Var1','num')

a$Var1= factor(a$Var1,levels = factor(order(-(as.numeric(a$Var1)))))
a$group <-if_else(a$Var1 ==1,'Group1','NA')
a$group <-if_else(a$Var1 ==2,'Group2',a$group)
a$group <-if_else(a$Var1 ==3,'Group3',a$group)
a$group <-if_else(a$Var1 == 4 | a$Var1 == 5,'Group4',a$group)
a$group <-if_else(a$Var1 == 6 | a$Var1 == 7 | a$Var1 == 8
                       | a$Var1 == 9 | a$Var1 == 10,'Group5',a$group)
a$group <-if_else(a$group =='NA','Group6',a$group)

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)

library(ggplot2)
p <- ggplot(data = a,aes(x = Var1,y= num)) +
  geom_col(aes(fill = group),position=position_dodge(1),width = 1,color = 'white') +
  coord_flip(ylim = c(0,140000)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=num,y=num + 3500), color="black", size=3) +
  cowplot::theme_half_open()+
  scale_fill_manual(values = b) +
  ylab("Number of eG4s ") + 
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.ticks.x = element_blank(),
        legend.position = "none")

#3500
Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/number_of_G4.pdf",
                  width = 8,height = 5,devices = NULL,res = 300,show = F) 

#4500
Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/number_of_G4.2.pdf",
                  width = 6.8,height = 6.8,devices = NULL,res = 300,show = F)
# -------------------------------------------------------------------------
hg_pqs$group <-if_else(hg_pqs$sum ==0,'PQS','NA')
hg_pqs$group <-if_else(hg_pqs$sum ==1,'Group1',hg_pqs$group)
hg_pqs$group <-if_else(hg_pqs$sum ==2,'Group2',hg_pqs$group)
hg_pqs$group <-if_else(hg_pqs$sum ==3,'Group3',hg_pqs$group)
hg_pqs$group <-if_else(hg_pqs$sum == 4 | hg_pqs$sum == 5,'Group4',hg_pqs$group)
hg_pqs$group <-if_else(hg_pqs$sum == 6 | hg_pqs$sum == 7 | hg_pqs$sum == 8
                       | hg_pqs$sum == 9 | hg_pqs$sum == 10,'Group5',hg_pqs$group)
hg_pqs$group <-if_else(hg_pqs$group =='NA','Group6',hg_pqs$group)

#write.table(hg_pqs,file = '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.group.result',
sep = '\t',row.names = F,quote = F)

# 4.number_of_G4.group-------------------------------------------------------------------------
hg_pqs <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.group.result') %>% data.frame()
b <- table(hg_pqs$group) %>% data.frame()
colnames(b) <- c('Var1','num')
b <- b[-7,]

b$Var1 %<>% factor(.,levels = factor(c(paste0('Group',1:6))))

a <- colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
#colors =  c('#fddbc7','#f4a582','#DB8872','#962E20','#b2182b','#67001f')
#a1 <- c("#EEDDD3","#D7AF9E","#AE7369",'#d6604d',"#7F323C","#400D1C")

p <- ggplot(data = b,aes(x = Var1,y= (num/100000),fill = Var1)) +
  geom_bar(stat = 'identity',position=position_dodge(0.7),width = 0.7,color = 'white') +
  coord_cartesian(ylim = c(0,1.3)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = a) +
  geom_text(aes(label=num,vjust = -0.5), color="black", size=4) +
  cowplot::theme_half_open()+
  ylab(bquote(Number~of~G4~(x10^5))) + 
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14),
        legend.position = "none")

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/number_of_G4.group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# 5.structure_stablity.7groups  --------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.group.result') %>% data.frame()
df$group %<>% gsub('PQS','Non-eG4',.)
#计算中位数
medians <- aggregate(score ~  group, df, median)
a <- table(df$group) %>% data.frame()
colnames(a)[1] <- 'group'
colnames(a)[2] <- 'number'
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
my_comparisons = list(c("Non-eG4","Group1"),c("Group1","Group2"),
                      c("Group2","Group3"),c("Group3","Group4"),
                      c("Group4","Group5"),c("Group5","Group6"))

library(ggpubr)
df$group %<>% factor(.,levels = c("Non-eG4",paste0('Group',1:6)))

p <- ggplot(data = df,aes(x = group,y = score,fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  #geom_text(data = medians, aes(label = score, y = 40),size = 4,col = 'red')+
  #geom_text(data = a,aes(label = as.character(number),y=45),size = 4,col = 'black')+
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(105,115),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("Structural stability")+
  coord_cartesian(ylim=c(40,150))+
  geom_hline(yintercept = medians[medians$group =='Non-eG4',2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/structure_stablity.7groups.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 
# density.G4_length-------------------------------------------------------------------------
G4 <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed') %>% data.frame()
G4$length <- G4$V3-G4$V2
n <- vector()
n[2] <- mean(G4$length)
n[1] <- nrow(G4)
p <- ggplot(G4,aes(length))+geom_bar(aes(y = (..count..)/sum(..count..)*100),width = 0.8,fill="#D7AF9E")+
  xlab("G4 length (bp)") + 
  ylab("Density (%)")+
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text= element_text(size=14,colour = "black"),
        axis.title = element_text(size=16))+
  scale_y_continuous(expand = c(0,0),limits = c(0,5))+
  annotate(geom="text", x=38, y=4.3, size=4,
           label=paste0("Total G4: ",n[1],"\n","Mean length: ",ceiling(n[2])," bp"))

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/density.G4_length.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)    
# length.PQS_G4-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
G4 <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed') %>% data.frame()
pqs <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed') %>% data.frame()

pqs$length <- pqs$V3-pqs$V2
colnames(pqs)[5] <- 'score'
pqs$type <- 'PQS'

G4$length <- G4$V3-G4$V2
colnames(G4)[5] <- 'score'
G4$type <- 'G4'

df <- rbind(pqs,G4)

#计算中位数
medians <- aggregate(length ~ type , df, median)

library(ggplot2)
library(ggpubr)
p <- ggplot(data = df,aes(x = type,y = length,fill = type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white",alpha = 0.8) +
  scale_x_discrete(limits=c("PQS","G4")) +
  scale_fill_manual(values = c('#b2182b','#4d4d4d')) +
  stat_compare_means(comparisons = list(c("PQS","G4")),label.y = 55,
                     aes(label = paste0("p = ", ..p.format..)),
                     tip.length = 0,size=5)+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("Length (bp)")+
  coord_cartesian(ylim=c(10,60))+
  geom_hline(yintercept = medians[2,2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/length.PQS_G4.pdf",
                  width = 4,height = 3,devices = NULL,res = 300,show = F)   
# feature.frequency.facet_wrap-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(magrittr)
path <- "/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

peakAnno <- list()
group <- list()
for (i in c(1:length(files))){
  peakAnno[[i]] <- annotatePeak(filepath[[i]], tssRegion=c(-1500, 1500),
                                TxDb=txdb, annoDb="org.Hs.eg.db")  
}
for (i in c(1:length(files))){
  group[[i]] <- peakAnno[[i]]@annoStat
}

a <- lapply(names(filepath), function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
a[7] <- 'non-eG4'
library(Hmisc)
a <- capitalize(a)
names(group) <- a
for (i in c(1:length(files))){
  group[[i]]$group <- a[i]
}

df <- do.call('rbind',group)
df$group %<>% factor(.,c('Non-eG4',paste0('Group',c(1:6))))

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

p <- ggplot(data = df,aes(x = Feature,y= Frequency)) +
  geom_col(aes(fill = group),position = 'stack',width = 1,
           color ="white") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()  +
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust=1,hjust= 1,size=12),
        axis.text.y = element_text(size=10),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  facet_wrap(.~group,scales= "free_y",ncol = 2)+
  ylab("Frequency (%)")+
  geom_text(aes(label=str_c(round(Frequency,2)),hjust = -0.1),
            color="black", size=3) +
  coord_flip(ylim = c(0,100))

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/feature.frequency.facet_wrap.pdf",
                  width = 6.5,height = 7,devices = NULL,res = 300,show = F) 
# read_count_frequency.genomic_region-------------------------------------------------------------------------
rm(list = ls());gc()
path <- "/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

#加载包
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

peak=list(`Non-eG4` = filepath[[7]],Group1= filepath[[1]],Group2 = filepath[[2]],
          Group3 = filepath[[3]],Group4 = filepath[[4]],Group5 = filepath[[5]],
          Group6 = filepath[[6]])

promoter <- getPromoters(TxDb=txdb, upstream=1500, downstream=1500)
tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter)

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

p <- plotAvgProf(tagMatrixList, xlim=c(-1500, 1500),
                 xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") +
  cowplot::theme_half_open()+
  theme(axis.title = element_text(size= 16),
        axis.text = element_text(size= 14),
        legend.position = "bottom") +labs(color = 'Group')+
  scale_fill_manual(values = color,aesthetics = "color")
                    
Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/read_count_frequency.genomic_region.pdf",
                  width = 7.5,height = 5,devices = NULL,res = 300,show = F) 
# GC_content.G4_length-------------------------------------------------------------------------
rm(list = ls());gc()
library(data.table)
library(magrittr)
G4_length <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4_length')
gc <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/1k_gc.bed')


df <- G4_length[,c(1:3,10)]
df$V5 <- paste(df$V1,df$V2,df$V3,sep = '_')
df2 <- tapply(df$V10,df$V5,sum) %>% data.frame()
colnames(df2)[1] <- 'length'
df2 <- rownames_to_column(df2 )

colnames(gc) <- c('V1','V2','V3','GC')
gc$V5 <- paste(gc$V1,gc$V2,gc$V3,sep = '_')

df2$GC <- gc[match(df2$rowname,gc$V5),4] %>% unlist()

p <- ggplot(data = df2,aes(x = length,y = GC)) +
  geom_point(shape = 21,size = 0.5,color = '#AE7369',fill = 'white', alpha = 0.2) +
  cowplot::theme_half_open()+
  xlab("G4 length (bp)")+
  ylab("GC content (%)")+
  annotate("text",x = 200,y = 0.1,label = "r = 0.3; pvalue = 2.2e-16",size = 5) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "none")

cor.test(df2$length,df2$GC)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/GC_content.G4_length.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)     
# GO-------------------------------------------------------------------------
rm(list = ls());gc()
library(clusterProfiler)
df <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/go.txt') %>% data.frame()

ann_col<-c("PQS",paste0('Group',1:6))
names(ann_col) <- ann_col
ann_col<-data.frame(group=ann_col)

ann_colors = list(
  group = c(PQS="#4d4d4d",Group1="#FDDBC7",Group2="#DFAFA5",
            Group3="#C18383",Group4="#A35762",Group5="#842B40",Group6 = '#67001F'))

p <- pheatmap::pheatmap(df,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        annotation_col =ann_col,
                        cellwidth=20,
                        cellheight = 10,
                        show_colnames = FALSE,
                        annotation_names_col = FALSE,
                        annotation_names_row = FALSE,
                        color = colorRampPalette(colors = c("white","#08306b"))(100),
                        annotation_colors = ann_colors,
                        legend_breaks =c(0,10,20,30))

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/GO.pdf",
                  width = 8,height = 3,devices = NULL,res = 300,show = F) 
# proportion_in_state.genome_region.group-------------------------------------------------------------------------
data <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/chromHMM.Type.ratio.txt') %>% data.frame()

color <- c("#FF0000","#FF4500","#32CD32","#008000","#006400","#C2E105","#FFFF00","#66CDAA",
           "#8A91D0","#CD5C5C","#E9967A","#BDB76B","#808080","#C0C0C0","#e5e5e5")
data$Group.1 %<>% factor(.,levels =  c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx", "5_TxWk",       
                                       "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv",
                                       "11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies"))
data$Group.2 %<>%  gsub('G4','eG4',.) %>% gsub('PQS','Non-eG4',.)
data$Group.2 %<>% factor(.,levels = c('Non-eG4',"eG4","Genome","Promoter","5'UTR","CDS","3'UTR",
                                      "Exon","Intron",paste0('Group',1:6)))
library(ggplot2)

p <- ggplot(data ,aes(x = Group.2, y = ratio,fill = Group.1)) +
  geom_bar(stat = 'identity',position = 'stack',width= 0.5) +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  ylab('Proportion in state')+
  theme(axis.title.x =element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size=14,angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size=14),
        legend.position = 'bottom')+
  geom_vline(xintercept = 9.5,linetype = 2,size=1.0)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/proportion_in_state.genome_region.group.pdf",
                  width = 8.2,height = 6,devices = NULL,res = 300,show = F) 
# state.proportion-------------------------------------------------------------------------
data <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/state.overlap.ratio.txt') %>% data.frame()

color <- rev(c('black',"#FF0000","#FF4500","#32CD32","#008000","#006400","#C2E105","#FFFF00","#66CDAA",
               "#8A91D0","#CD5C5C","#E9967A","#BDB76B","#808080","#C0C0C0","#e5e5e5"))

data$state %<>% factor(.,levels =  rev(c('Bases',"1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx", "5_TxWk",       
                                         "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv",
                                         "11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")))

p <- ggplot(data ,aes(x = state, y = (ratio*100),fill = state)) +
  geom_bar(stat = 'identity',width=0.7) +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  ylab('Proportion of state overlapping eG4s (%)')+
  theme(axis.title.x =element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14),
        legend.position = 'none') +
  coord_flip(ylim = c(0,10))+
  geom_text(data = data, aes(label = round(ratio*100, 2), y = ratio*100 + 0.5), size = 4) +
  geom_vline(xintercept = 15.5,linetype=2,size = 0.8)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/state.proportion.pdf",
                  width = 7,height = 4.5,devices = NULL,res = 300,show = F) 
#Proportion_in_each_G4_group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/Proportion_in_each_G4_group.txt') %>% data.frame()
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)

df$Group.1 %<>% factor(.,levels = rev(c('Bases',"1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx", "5_TxWk",       
                                           "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv",
                                           "11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")))

p <- ggplot(df ,aes(x = Group.1, y = ratio,fill = Group.2)) +
  geom_bar(stat = 'identity',width=0.8,color = 'white') + 
  scale_fill_manual(values = b) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  ylab('Proportion in each eG4 group')+
  theme(axis.title =element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'bottom') +
  coord_flip()+
  geom_vline(xintercept = 15.5,linetype=2)
Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/Proportion_in_each_G4_group.pdf",
                  width = 5,height = 4.5,devices = NULL,res = 300,show = F) 
# DHS.genome_region.group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/DHS.genome_region.group.txt') %>% data.frame()
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
library(scales)
mypal <- pal_npg("nrc", alpha=1)(9)
mypal <- mypal[2:7]
col <- c('#4d4d4d','#b2182b',mypal,b)

df$type %<>%  gsub('G4','eG4',.) %>% gsub('PQS','Non-eG4',.)
df$type %<>% factor(.,levels = c('Non-eG4',"eG4","Promoter","5'UTR","CDS","3'UTR",
                                      "Exon","Intron",paste0('Group',1:6)))
my_comparisons = list(c("Non-eG4","eG4"),c("eG4","Promoter"),
                      c("Promoter","5'UTR"),c("5'UTR","CDS"),
                      c("CDS","3'UTR"),c("3'UTR",'Exon'),
                      c('Exon','Intron'),c('Group1','Group2'),
                      c('Group2','Group3'),c('Group3','Group4'),
                      c('Group4','Group5'),c('Group5','Group6'))

library(ggsci)
library(ggpubr)
p <- ggplot(data = df,aes(x = type,y = ratio,fill = type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = col) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(1,1.1),7)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size= 3)+
  cowplot::theme_half_open() +
  theme(axis.text.y = element_text(size= 14),
        axis.text.x = element_text(size= 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab(paste0("Proportion of genomic region",'\n',"overlapping DHS"))+
  geom_vline(xintercept = 8.5,linetype=2,size = 0.6)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/DHS.genome_region.group.pdf",
                  width = 7,height = 5,devices = NULL,res = 300,show = F) 
# H3K27ac.genome_region.group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/H3K27ac.genome_region.group.txt') %>% data.frame()
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
library(scales)
library(ggsci)
library(ggpubr)
mypal <- pal_npg("nrc", alpha=1)(9)
mypal <- mypal[2:7]
col <- c('#4d4d4d','#b2182b',mypal,b)

df$type %<>%  gsub('G4','eG4',.) %>% gsub('PQS','Non-eG4',.)
my_comparisons = list(c("Non-eG4","eG4"),c("eG4","Promoter"),
                      c("Promoter","5'UTR"),c("5'UTR","CDS"),
                      c("CDS","3'UTR"),c("3'UTR",'Exon'),
                      c('Exon','Intron'),c('Group1','Group2'),
                      c('Group2','Group3'),c('Group3','Group4'),
                      c('Group4','Group5'),c('Group5','Group6'))

df$type %<>% factor(.,levels = c("Non-eG4","eG4","Promoter","5'UTR","CDS","3'UTR",
                                 "Exon","Intron",paste0('Group',1:6)))

p <- ggplot(data = df,aes(x = type,y = ratio,fill = type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = col) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(1,1.1),7)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size= 3)+
  cowplot::theme_half_open() +
  theme(axis.text.y = element_text(size= 14),
        axis.text.x = element_text(size= 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab(paste0("Proportion of genomic region",'\n',"overlapping H3K27ac")) +
  geom_vline(xintercept = 8.5,linetype=2,size = 0.6)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/H3K27ac.genome_region.group.pdf",
                  width = 7,height = 5,devices = NULL,res = 300,show = F) 
# number_of_chromHMM_states-------------------------------------------------------------------------
d <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/number_of_chromHMM_states.txt') %>% data.frame()
d$state %<>% factor(.,levels = c(1:15))
d$Group %<>% factor(.,levels = c('PQS',paste0('Group',1:6))) 

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
col <- c('#4d4d4d',b)

ggplot(d,aes(x = state,y = ratio,color = Group)) +
  geom_line(aes(group = Group),size = 1) +
  geom_point(shape = 21,size = 3,fill = "white") + #21空心点
  scale_color_manual(values = col) +
  cowplot::theme_half_open() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        #axis.ticks.x = element_blank(),
        legend.position = c(0.8,0.7)) +
  xlab("Number of chromHMM states")+
  ylab("Fraction of G4s (%)")
# conservation-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata2/conservation_score/conservation_information') %>% data.frame()
df$group %<>% gsub('PQS','non-eG4',.)
df$group %<>% factor(.,levels = c("non-eG4",paste0('Group',1:6)))
#计算中位数
medians <- aggregate(phastCons ~  group, df, median)

my_comparisons = list(c("non-eG4","Group1"),c("Group1","Group2"),
                      c("Group2","Group3"),c("Group3","Group4"),
                      c("Group4","Group5"),c("Group5","Group6"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

p <- ggplot(data = df,aes(x = group,y = phastCons,fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white",size = 0.6) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(0.9,1),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size = 4)+
  scale_fill_manual(values = color ) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  coord_cartesian(ylim = c(-0.1,1.2))+
  ylab("PhastCons Score")+
  geom_hline(yintercept = medians[medians$group =='non-eG4',2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/phastCons_score.group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

medians2 <- aggregate(phyloP ~  group, df, median)
p2 <- ggplot(data = df,aes(x = group,y = phyloP,fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white",size = 0.6) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(2,1.7,2,1.7,2,1.7),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color ) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  coord_cartesian(ylim = c(-2,3))+
  ylab("PhyloP Score")+
  geom_hline(yintercept = medians2[medians2$group =='non-eG4',2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/phyloP_score.group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

df2 <- df[df$group != 'PQS',]
PQS <- df[df$group == 'PQS',]
n = c((quantile((PQS$phastCons),0.95) %>% round(.,3)),(quantile((PQS$phyloP),0.95) %>% round(.,3)))
res <- mutate(df2,direction = if_else(abs(phastCons) > n[1] & abs(phyloP) > n[2],'C','NS'))

p3 <- ggplot(data = res,aes(x = phastCons,y = phyloP)) +
  geom_point(shape = 21,size = 0.5,aes(fill = direction),alpha = 0.08) +
  cowplot::theme_half_open() +
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  xlab("PhastCons Score")+
  ylab("PhyloP Score")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))+
  geom_hline(yintercept = n[2],
             linetype = 2,size = 0.7) +
  geom_vline(xintercept = n[1],
             linetype = 2,size = 0.7) +
  scale_fill_manual(values = c('#b2182b','#4d4d4d'))+
  coord_cartesian(ylim = c(-2.5,7.5))

Ipaper::write_fig(p3,file = "/home/yulix/G4_analysis/result/figure2/phastCons_score.phyloP_score.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 


a <- c(paste0('Group',1:6))
prop= vector()
for (i in a){
  prop[i] <- sum(df2$group == i & abs(df2$phastCons) > n[1] & abs(df2$phyloP) > n[2])/sum(df2$group ==i)
}
data <- data.frame(prop)
data <- rownames_to_column(data)
colnames(data)[1] <- 'group'

data$group %<>% factor(.,levels = c(paste0('Group',6:1)))
color2 <- colorRampPalette(colors = c('#67001f', '#fddbc7'))(6)

p4 <- ggplot(data ,aes(x = group, y = prop,fill = group)) +
  geom_bar(stat = 'identity',width=0.8,color = 'white') +
  geom_text(color = 'black',
            aes(x = group,y=prop,hjust = 1,
                label = percent(prop,digits = 3)),size=5) +
  scale_fill_manual(values=color2) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  ylab(' Proportion of conservation G4')+
  theme(axis.title.y =element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = 'none') +
  coord_flip(ylim =c(0,0.35))

Ipaper::write_fig(p4,file = "/home/yulix/G4_analysis/result/figure2/group.proportion_of_conservation.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# RNA-seq-------------------------------------------------------------------------
# fraction.G4_group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/fraction.G4_group.txt') %>% data.frame()
df$Class %<>% gsub('G4-depleted genes','eG4-depleted genes',.) %>% gsub('G4-containing genes','eG4-containing genes',.)

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
df$Class %<>% factor(.,levels = c('eG4-depleted genes','eG4-containing genes',c(paste0('Group',1:6))))  

p <- ggplot(data = df,aes(x = Class, y = Fraction,fill = Class)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab('Fraction of expressed genes') +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,1))+
  geom_text(data = df, aes(label = round(Fraction,3)),vjust = -1, size = 5)+
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/fraction.G4_group.pdf",
                  width = 7.5,height = 6.5,devices = NULL,res = 300,show = F)        

# expression_level.G4_group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/expression_level.G4_group.txt') %>% data.frame()
df$class %<>% gsub('G4-depleted genes','eG4-depleted genes',.) %>% gsub('G4-containing genes','eG4-containing genes',.)

#medians.exp <- aggregate(expression ~ class,df,median)

df$class %<>% factor(.,levels = c('eG4-depleted genes','eG4-containing genes',c(paste0('Group',1:6))))  
my_comparisons <- list(c('eG4-depleted genes','eG4-containing genes'),c('Group1','Group2'),c('Group2','Group3'),
                       c('Group3','Group4'),c('Group4','Group5'),c('Group5','Group6'))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
library(ggpubr)
p <- ggplot(data = df,aes(x = class,y = log10(expression+1),fill = class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white")+
  coord_cartesian(ylim = c(0,6)) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(5,4.2),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  cowplot::theme_half_open()+
  scale_fill_manual(values = color) +
  ylab(Expression~level~(log[10](TPM))) +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position = "none" ) +
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/expression_level.G4_group.pdf",
                  width = 7.5,height = 6.5,devices = NULL,res = 300,show = F)        

# promoter_fraction.G4_group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/promoter_fraction.G4_group.txt') %>% data.frame()

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
df$Class %<>% factor(.,levels = c('G4-depleted genes','G4-containing genes',c(paste0('Group',1:6))))  

p <- ggplot(data = df,aes(x = Class, y = Fraction,fill = Class)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab('Fraction of expressed genes') +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,1))+
  geom_text(data = df, aes(label = round(Fraction,3)),vjust = -1, size = 5)+
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/promoter_fraction.G4_group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)    

# promoter_expression_level.G4_group-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/promoter_expression_level.G4_group.txt') %>% data.frame()

medians.exp <- aggregate(exp ~ group,df,median)

df$group %<>% factor(.,levels = c('G4-depleted genes','G4-containing genes',c(paste0('Group',1:6))))  
my_comparisons <- list(c('G4-depleted genes','G4-containing genes'),c('Group1','Group2'),c('Group2','Group3'),
                       c('Group3','Group4'),c('Group4','Group5'),c('Group5','Group6'))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
library(ggpubr)
p <- ggplot(data = df,aes(x = group,y = log10(exp+1),fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white")+
  coord_cartesian(ylim = c(0,6)) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(5,4.2),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  cowplot::theme_half_open()+
  scale_fill_manual(values = color) +
  ylab(Expression~level~(log[10](TPM))) +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position = "none" ) +
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/promoter_expression_level.G4_group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# enhancer-------------------------------------------------------------------------
df = fread('/home/yulix/G4_analysis/result/savedata/epigenetics/enhancer.txt')%>% data.frame() 
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
df$group %<>% factor(.,levels = c('PQS','G4',c(paste0('Group',1:6))))  

ggplot(data = df,aes(x = group, y = ratio,fill = group)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab('Fraction of expressed genes') +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,1))+
  geom_text(data = df, aes(label = round(ratio,2)),vjust = -1, size = 5)+
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)+facet_wrap(.~class)

# overlap_length.strand-------------------------------------------------------------------------
df <- data.table::fread('/home/yulix/G4_analysis/result/savedata/gene_region_intersect_G4/overlap_length.strand.txt') %>% data.frame()
df$type %<>% factor(.,levels = c('PQS','G4',paste0('Group',1:6)))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)

ggplot(df ,aes(x = class,y = mean,fill = strand)) +
  geom_bar(stat="identity",width= 0.75,position= position_dodge(0.75),color = 'white') +
  #scale_fill_manual(values = color) +
  scale_fill_manual(values = c('#4d4d4d','#b2182b')) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  ylab('')+
  theme(axis.title.y =element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = 'top') +
  facet_wrap(.~type)
#coord_cartesian(ylim = c(0,12))

ggplot(df ,aes(x = class,y = mean,fill = type)) +
  geom_bar(stat="identity",width= 0.75,position= position_dodge(0.75),color = 'white') +
  scale_fill_manual(values = color) +
  #scale_fill_manual(values = c('#4d4d4d','#b2182b')) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  ylab('')+
  theme(axis.title.y =element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.position = 'top') +
  facet_wrap(.~strand)

# overlap_num.strand-------------------------------------------------------------------------
df <- data.table::fread('/home/yulix/G4_analysis/result/savedata/gene_region_intersect_G4/overlap_num.strand.txt') %>% data.frame()
df$type <- if_else(df$type == 'G4','eG4',df$type)
df1 <- filter(df,type %in% c('Non-eG4','eG4'))
df2 <- filter(df,type %in% c(paste0('Group',1:6)))

df1$type %<>% factor(.,levels = c('Non-eG4','eG4'))


p <- ggplot(df1 ,aes(x = class,y = mean,fill = strand)) +
  geom_bar(stat="identity",width= 0.75,position= position_dodge(0.75),color = 'white') +
  scale_fill_manual(values = c('#4d4d4d','#b2182b')) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  ylab(paste0('Density of eG4 and non-eG4','\n','(number/10^4 base pair)'))+
  theme(axis.title.y =element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = 'top') +
  facet_wrap(.~type,scales = 'free')

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/overlap_num.strand.non-eG4_eG4.pdf",
                  width = 7.5,height = 5,devices = NULL,res = 300,show = F) 

p2 <- ggplot(df2 ,aes(x = class,y = mean,fill = strand)) +
  geom_bar(stat="identity",width= 0.75,position= position_dodge(0.75),color = 'white') +
  scale_fill_manual(values = c('#4d4d4d','#b2182b')) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  ylab(paste0('Density of eG4 and non-eG4','\n','(number/10^4 base pair)'))+
  theme(axis.title.y =element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = 'top') +
  facet_wrap(.~type,scales = 'free')

Ipaper::write_fig(p2,file = "/home/yulix/G4_analysis/result/figure2/overlap_num.strand.eG4_group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# fold_enrichment.PQS_G4-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/fold_enrichment.PQS_G4.txt') %>% data.frame()

df$group %<>% factor(.,levels = c('PQS','G4',paste0('Group',1:6)))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)

p <- ggplot(data = df,aes(x = group , y = FoldEnrichment)) +
  geom_bar(stat="identity",width=0.8,position='dodge',aes(fill = group),color = 'white') +
  scale_fill_manual(values = color )+
  geom_text(data = df, aes(label = round(FoldEnrichment,3)),vjust = -1,size = 4.5) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Fold enrichment\n(overlapped with TE)")+
  cowplot::theme_half_open() +
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = "none") +
  geom_hline(yintercept = 1,linetype = 2,size = 1)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/fold_enrichment.PQS_G4.pdf",
                  width = 7,height = 4,devices = NULL,res = 300,show = F) 
# fold_enrichment.desease-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/fold_enrichment.desease.txt') %>% data.frame()

df$class %<>% factor(.,levels = c('Lethal','Lethal_disease','Disease','OtherGenes','Oncogene','TSG'))
df$group %<>% gsub('G4','eG4',.) %>% gsub('PQS','Non-eG4',.)
df$group %<>% factor(.,levels = c('Non-eG4','eG4',paste0('Group',1:6)))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)

p <- ggplot(df,aes(x= class,y = Fold.enrichment,fill = group)) +
  geom_bar(stat = "identity",position=position_dodge(0.8),width=0.8,color = 'white') +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()+
  scale_fill_manual(values = color) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = c(0.3,0.9),legend.direction = 'horizontal') +
  coord_cartesian(ylim = c(0,2.5)) +
  ylab("Fold enrichment")+
  geom_vline(xintercept = 4.5,linetype = 1) +
  geom_hline(yintercept = 1,linetype = 2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/fold_enrichment.desease.pdf",
                  width = 7,height = 4,devices = NULL,res = 300,show = F)   

df[df$group == 'G4',4]/df[df$group == 'PQS',4]
factor(df$class)

# promoter_fold_enrichment.desease-------------------------------------------------------------------------

df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/promoter_fold_enrichment.desease.txt') %>% data.frame()
df$class %<>% factor(.,levels = c('Lethal','Lethal_disease','Disease','OtherGenes','Oncogene','TSG'))
df$group %<>% factor(.,levels = c('PQS','G4',paste0('Group',1:6)))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)

ggplot(df,aes(x= class,y = Fold.enrichment,fill = group)) +
  geom_bar(stat = "identity",position=position_dodge(0.8),width=0.8,color = 'white') +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()+
  scale_fill_manual(values = color) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = 'none') +
  #coord_cartesian(ylim = c(0,2.5)) +
  ylab("Fold enrichment")+
  geom_vline(xintercept = 4.5,linetype = 1) +
  geom_hline(yintercept = 1,linetype = 2)

df[df$group == 'G4',4]/df[df$group == 'PQS',4]
factor(df$class)

# density.G4_desease-------------------------------------------------------------------------
G4 <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/intersect_obs/G4.bed') %>% data.frame()
G4 <- G4[G4$V9 != 0,]
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_tsg_oncogene.txt') %>% data.frame()
G4$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist() 
df$num <- G4[match(df$gene_id,G4$V4),9]
df$class %<>% gsub('Tsg','TSG',.)

df$class %<>% factor(.,levels = c('Lethal','Lethal_disease','Disease','OtherGenes','Oncogene','TSG'))

sum(is.na(df))
df <- na.omit(df)

a = aggregate(df$num,by = list(df$class),mean)
colnames(a) <- c('class','mean')

lab = c(lapply(a$mean,function(x)paste0("Mean Load","\n",round(x,0)," G4s")) %>% unlist())

ggplot(df,aes(num))+
  geom_bar(aes(y = (..count..)/sum(..count..)*100),width = 1,fill="#D7AF9E")+
  xlab("Number of G4") + 
  ylab("% Genes")+
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text= element_text(size=14,colour = "black"),
        axis.title = element_text(size=16))+
  scale_y_continuous(expand = c(0,0),limits = c(0,5))+
  coord_cartesian(xlim = c(0,100)) +
  geom_text(aes(x,y,label = lab),data = data.frame(x=rep(50,6),y = rep(3,6),lab = lab,class = levels(df$class)),size = 4,vjust = 1)+
  facet_wrap(.~class)

df[,4:6] <- G4[match(df$gene_id,G4$V4),1:3]
df <- df[,c(4:6,1:3)] 
b <- unique(df$class) %>% as.vector()

list <- list()
for (i in 1:length(b)){
  list[[i]] <- df[df$class == b[i],]
}
names(list) <- b

#批量导出结果
outpath <- '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/class'
out_filename <- sapply(names(list),function(x){
  paste(x,'.tsv',sep='')})
outfilepath <- sapply(out_filename,function(x){paste(outpath,x,sep='/')})
for (i in 1:length(list)){
  write.table(list[[i]],file = outfilepath[[i]],sep = '\t',row.names = FALSE,col.names = F,quote=FALSE)
}


# read_count_frequency.desease-------------------------------------------------------------------------

rm(list = ls());gc();rm(list = ls())
path <- "/home/yulix/G4_analysis/result/savedata/G4_TE/disease/class"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

#加载包
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

peak=list(Lethal = filepath[[3]],Lethal_disease = filepath[[2]],Disease = filepath[[1]],
          Oncogene = filepath[[4]],TSG = filepath[[6]])

promoter <- getPromoters(TxDb=txdb, upstream=1500, downstream=1500)
tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter)

b = colorRampPalette(colors = c('white', '#b2182b'))(4)
color <- c(b[4:2],'#4393c3','#92c5de')

p <- plotAvgProf(tagMatrixList, xlim=c(-1500, 1500),
                 xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") +
  cowplot::theme_half_open()+
  theme(axis.title = element_text(size= 16),
        axis.text = element_text(size= 14),
        legend.position = "bottom") +labs(color = 'Class')+
  scale_fill_manual(values = color,aesthetics = "color")

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/read_count_frequency.desease.pdf",
                  width = 7,height = 5,devices = NULL,res = 300,show = F) 