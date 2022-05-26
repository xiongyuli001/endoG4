
# feature.frequency.facet_wrap-------------------------------------------------------------------------
pqs <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.group.result') %>% data.frame()
group1 <- pqs[pqs$group == 'Group1',1:6]
group2 <- pqs[pqs$group == 'Group2',1:6]
group3 <- pqs[pqs$group == 'Group3',1:6]
group4 <- pqs[pqs$group == 'Group4',1:6]
group5 <- pqs[pqs$group == 'Group5',1:6]
group6 <- pqs[pqs$group == 'Group6',1:6]

list = list(group1 = group1,group2 = group2,group3 = group3,
         group4 = group4,group5 = group5,group6 = group6)

#批量导出
outpath <- '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group'
out_filename <- sapply(names(list),function(x){
  paste(x,'.bed',sep='')})
outfilepath <- sapply(out_filename,function(x){paste(outpath,x,sep='/')})
for (i in 1:length(list)){
  write.table(list[[i]],file = outfilepath[[i]],sep = '\t',row.names = FALSE,col.names = F,quote=FALSE)
}
# GO-------------------------------------------------------------------------
rm(list = ls());gc()
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
path <- "/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})
peak=list(PQS = filepath[[7]],Group1= filepath[[1]],Group2 = filepath[[2]],
          Group3 = filepath[[3]],Group4 = filepath[[4]],Group5 = filepath[[5]],
          Group6 = filepath[[6]])
peakAnnoList = list()
for (i in 1:length(peak)) {
  peakAnnoList[[i]] <- lapply(peak[[i]], annotatePeak, TxDb=txdb,
                              tssRegion=c(-1500, 1500), verbose=FALSE)
}         
gene=list()
for (i in 1:length(peak)) {
  gene[[i]] = lapply(peakAnnoList[[i]], function(i) as.data.frame(i)$geneId)
}
list <- list()
for (i in 1:length(peak)) {
  list[[i]] <- data.frame(t(data.frame(gene[[i]])))
  list[[i]] <- unlist(as.character(list[[i]]))
}
names(list) <- c("PQS",paste0('Group',1:6))

erich.go.BP <- list()
for (i in c(1:7)) {
  erich.go.BP[[i]] = enrichGO(gene = list[[i]],
                              OrgDb = org.Hs.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1)
}
#批量导出GO BP结果
outpath <- '/home/yulix/G4_analysis/result/savedata2/GO_BP'
out_filename <- sapply(names(list),function(x){
  paste(x,'.tsv',sep='')})
outfilepath <- sapply(out_filename,function(x){paste(outpath,x,sep='/')})
for (i in 1:length(erich.go.BP)){
  write.table(erich.go.BP[[i]],file = outfilepath[[i]],sep = '\t',row.names = FALSE,quote=FALSE)
}

a <- c("PQS",paste0('Group',1:6))
data = erich.go.BP

for (i in c(1:length(data))){
  data[[i]] <- arrange(data[[i]],p.adjust)
  data[[i]] <- data[[i]][1:5,]
  data[[i]] <- data[[i]][,c("Description","p.adjust")]
  data[[i]]$group <- a[i]
}

df <- do.call('rbind',data)
df <- dplyr::mutate(df,p.value = -log10(p.adjust)) 
df <- df[,c(1,3,4)]
max(df$p.value)
min(df$p.value)
df<-spread(
  data=df,   
  key=group,     
  value=p.value) 
df[is.na(df)]<- 0
library(Hmisc)
df$Description <- capitalize(df$Description)
df<-as.data.frame(column_to_rownames(df,var="Description"))
df <- df[,c(7,c(1:6))]

write.table(df,file = '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/go.txt',
            sep = '\t',col.names = T,row.names = T,quote = F)
# -------------------------------------------------------------------------
rm(list = ls());gc()
library(readr)
library(tidyverse)
gencode_v38lift37_annotation <- read_table2("/home/yulix/G4_analysis/data/ref/gencode.v38lift37.annotation.gtf")
a <- gencode_v38lift37_annotation[-c(1:4),]
b <- paste0('v',c(1:17))
colnames(a) <- b

gene <- a[a$v3 == 'gene',c(1,4,5,10,7)]
gene <- data.frame(gene[,1:4],'20',gene[,5])
gene$v10 <- lapply(gene$v10, function(x)gsub("\"","", x)) %>% unlist
gene$v10 <- lapply(gene$v10, function(x)gsub(";","", x)) %>% unlist
# promoter
gene_1 <- gene[gene$v7 == '+',]
gene_2 <- gene[gene$v7 == '-',]
promoter_1 <- gene_1[,c(1,2,4,6)]
promoter_2 <- gene_2[,c(1,3,4,6)]
promoter_1 <- data.frame(promoter_1$v1,promoter_1$v4-1500,
                         promoter_1$v4+1500,promoter_1$v10,promoter_1$v7)
colnames(promoter_1) <- c('chr','start','end','gene_id','strand')
promoter_1$start <- ifelse(promoter_1$start <0, 0, promoter_1$start)
promoter_2 <- data.frame(promoter_2$v1,promoter_2$v5-1500,
                         promoter_2$v5+1500,promoter_2$v10,promoter_2$v7)
colnames(promoter_2) <- c('chr','start','end','gene_id','strand')
promoter <- rbind(promoter_1,promoter_2)
write.table(promoter,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/promoter2.bed',sep = '\t',
            col.names = F,row.names = F,quote = F)

# -------------------------------------------------------------------------
rm(list = ls());gc()
library(data.table)
library(tidyverse)
path <- '/home/yulix/G4_analysis/result/savedata/epigenetics/genome_region'
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})
data <- list()
for (i in 1:length(files)){
  data[[i]] <- fread(filepath[[i]])
}
a = files
a <- lapply(a,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
a <- c("3'UTR","5'UTR","CDS","Exon","Genome", "Intron","Promoter") 

names(data) <- a
for (i in 1:length(files)){
  data[[i]] <- data[[i]][,1:3]
  data[[i]]$V4 <- a[i]
}
df <- do.call('rbind',data)

df2 <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.group.result') %>% data.frame()
df3 <- df2[,c(1:3,65)]
colnames(df3)[4] <- 'Type'
colnames(df) <- colnames(df3)

G4 <- df3[df3$Type != 'PQS',]
G4$Type <- 'G4'

df4 <- rbind(df,df3,G4)

write.table(df4,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/genoem_region.G4.pqs.group.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)


#bedtools intersect -a /home/yulix/G4_analysis/data/roadmap/input/chromHMM.hg19.merge.sort.bed -b genoem_region.G4.pqs.group.txt -wo > chromHMM_genome

#proportion_in_state.genome_region.group-----------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
df <- fread('~/G4_analysis/result/savedata/epigenetics/chromHMM_genome')
a = unique(df$V8)
b= unique(df$V4)

base <- vector()
for (i in 1:length(a)){
  base[i] = sum(df[df$V8 == a[i],]$V9)
}
base_num <- data.frame(a,base) 
data <- aggregate(x = df$V9,by = list(df$V4,df$V8),sum) %>% data.frame()
data$sum <- base_num[match(data$Group.2,base_num$a),2]
data$ratio <- data$x/data$sum

write.table(data,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/chromHMM.Type.ratio.txt',
            sep = '\t',col.names = T,row.names = F,quote = F)
# state.proportion-------------------------------------------------------------------------
rm(list = ls());gc()
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/chromHMM_genome')
chromHMM <- fread('/home/yulix/G4_analysis/data/roadmap/input/chromHMM.hg19.merge.sort.bed')

chrom_sum <- tapply(chromHMM$V3-chromHMM$V2, chromHMM$V4, sum) %>% data.frame()
chrom_sum <- rownames_to_column(chrom_sum)
colnames(chrom_sum) <- c('state','sum')

df <- df[df$V8 %in% c(paste0('Group',1:6)),]

overlap <- tapply(df$V9,df$V4, sum)%>% data.frame()
overlap <- rownames_to_column(overlap)
colnames(overlap) <- c('state','overlap')

data <- left_join(chrom_sum,overlap,by = 'state')
data$ratio <- data$overlap/data$sum
data[16,] <- c('Bases',0,0,0.004)
data$ratio %<>% as.numeric()

#G4 总长度占基因组比例
awk -F "\t" '{print $3-$2}' G4.bed|awk '{sum += $1};END {print sum}'
echo 12347722/3101804739|bc -l
.00398081860045800903

write.table(data,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/state.overlap.ratio.txt',
            sep = '\t',col.names = T,row.names = F,quote = F)
# Proportion_in_each_G4_group-------------------------------------------------------------------------
sed '1d' hg.pqs.group.result |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$NF}' |awk '{if($7!="PQS") print $0}' > 6groups.bed

rm(list = ls());gc()
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata/epigenetics/chromHMM.Type.ratio.txt') %>% data.frame()
df <- df[df$Group.2 %in% c(paste0('Group',1:6)),]

df$x %<>% as.numeric()
sum <- tapply(df$x, df$Group.1,sum) %>% data.frame()
sum <- rownames_to_column(sum)
colnames(sum) <- c('Group.1','x')
df$sum <- sum[match(df$Group.1,sum$Group.1),2]
df$ratio <- df$x/df$sum

G4 <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed') %>% data.frame()
G4$length <- G4$V3-G4$V2
a <- tapply(G4$length, G4$V7, sum) %>% data.frame()
a <- rownames_to_column(a)
colnames(a) <- c('Group.2','x')
sum(a$x) #12347722
a$sum <- 12347722
a$ratio <- a$x/a$sum
a$Group.1 <- 'Bases'

df2 <- rbind(a,df)

write.table(df2,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/Proportion_in_each_G4_group.txt',
            sep = '\t',col.names = T,row.names = F,quote = F)

# DHS.genome_region.group-------------------------------------------------------------------------
#7 个基因组区域,7 个G4 group和53 个DHS相交的数量所占的比例
for i in `ls /home/yulix/G4_analysis/result/savedata/epigenetics/genome_region/*bed` `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*bed`;do
for j in `ls /home/yulix/G4_analysis/data/roadmap/DHS/*DNase.macs2.narrowPeak`;do
a=`basename $i`
b=`basename $j`
bedtools intersect -a $i -b $j -c > ./dhs_num_ratio/${a%.bed}_${b%-DNase.macs2.narrowPeak}
done
done
#G4和53 个DHS相交的数量所占的比例
for j in `ls /home/yulix/G4_analysis/data/roadmap/DHS/*DNase.macs2.narrowPeak`;do
b=`basename $j`
bedtools intersect -a /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed -b $j -c > ./dhs_num_ratio/G4_${b%-DNase.macs2.narrowPeak}
done
#7 个基因组区域,7 个G4 group和98个h3k27ac相交的数量所占的比例
cd /home/yulix/G4_analysis/result/savedata/epigenetics
for i in `ls /home/yulix/G4_analysis/result/savedata/epigenetics/genome_region/*bed` `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*bed`;do
for j in `ls /home/yulix/G4_analysis/data/roadmap/H3K27ac/*narrowPeak`;do
a=`basename $i`
b=`basename $j`
bedtools intersect -a $i -b $j -c > ./h3k27ac_num_ratio/${a%.bed}_${b%-H3K27ac.narrowPeak}
done
done
#G4和98个h3k27ac相交的数量所占的比例
for j in `ls /home/yulix/G4_analysis/data/roadmap/H3K27ac/*narrowPeak`;do
b=`basename $j`
bedtools intersect -a /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed -b $j -c > ./h3k27ac_num_ratio/G4_${b%-H3K27ac.narrowPeak}
done
# DHS.genome_region.group-------------------------------------------------------------------------
rm(list = ls());gc()
library(data.table)
library(magrittr)
library(tidyverse)
path <- '/home/yulix/G4_analysis/result/savedata/epigenetics/dhs_num_ratio'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))

data <- list()
data <- lapply(filespath,function(x)fread(x))
names(data) <- files

b <- list()
for (i in 1:length(data)) {
  a =  data[[i]] %>% data.frame()
  b[[i]] <- a[,ncol(data[[i]])]
  b[[i]] <- if_else(b[[i]] >= 1,1,0)
}
names(b) <- files

c <- list()
for (i in 1:length(b)) {
  c[[i]] = sum(b[[i]] != 0)/length(b[[i]])
}
  
names(c) <- names(b)

df <- do.call('rbind',c) %>% data.frame()
df <- rownames_to_column(df)

df$V3 <- lapply(df$rowname,function(x)strsplit(x,'_E',fixed = T)[[1]][1]) %>% unlist()
colnames(df) <- c('V1','ratio','type')
df$type <- if_else(df$type == '5_UTR','5\'UTR',df$type)
df$type <- if_else(df$type == '3_UTR','3\'UTR',df$type)
df$type <- if_else(df$type == 'promoter2','promoter',df$type)
df$type <- if_else(df$type == 'pqs','PQS',df$type)
library(Hmisc)
df$type <- capitalize(df$type)
df <- df[df$type != 'Genome',]

write.table(df,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/DHS.genome_region.group.txt',
            sep = '\t',col.names = T,row.names = F,quote = F)
# H3K27ac.genome_region.group-------------------------------------------------------------------------
rm(list = ls());gc()
library(data.table)
library(magrittr)
library(tidyverse)
path <- '/home/yulix/G4_analysis/result/savedata/epigenetics/h3k27ac_num_ratio' 
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))

data <- list()
data <- lapply(filespath,function(x)fread(x))
names(data) <- files

b <- list()
for (i in 1:length(data)) {
  a =  data[[i]] %>% data.frame()
  b[[i]] <- a[,ncol(data[[i]])]
  b[[i]] <- if_else(b[[i]] >= 1,1,0)
}
names(b) <- files

c <- list()
for (i in 1:length(b)) {
  c[[i]] = sum(b[[i]] != 0)/length(b[[i]])
}

names(c) <- names(b)

df <- do.call('rbind',c) %>% data.frame()
df <- rownames_to_column(df)

df$V3 <- lapply(df$rowname,function(x)strsplit(x,'_E',fixed = T)[[1]][1]) %>% unlist()
colnames(df) <- c('V1','ratio','type')
df$type <- if_else(df$type == '5_UTR','5\'UTR',df$type)
df$type <- if_else(df$type == '3_UTR','3\'UTR',df$type)
df$type <- if_else(df$type == 'promoter2','promoter',df$type)
df$type <- if_else(df$type == 'pqs','PQS',df$type)
library(Hmisc)
df$type <- capitalize(df$type)
df <- df[df$type != 'Genome',]

write.table(df,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/H3K27ac.genome_region.group.txt',
            sep = '\t',col.names = T,row.names = F,quote = F)
# number_of_chromHMM_states-------------------------------------------------------------------------
#7 个G4 group和chromHMM_15_states相交的数量
cd ~/G4_analysis/result/savedata/epigenetics/chromHMM_15_states
for i in *;do
for j in /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*;do 
a=`basename $j`
bedtools intersect -a $j -b $i -c > ../7group_c/${i%bed}${a%.bed}
done
done

rm(list = ls());gc()
library(data.table)
library(tidyverse)
library(dplyr)
path <- '/home/yulix/G4_analysis/result/savedata/epigenetics/7group_c'
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})
df <- lapply(filepath,function(x)fread(x))  

a = lapply(names(df),function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
b = lapply(names(df),function(x)strsplit(x,'.',fixed = T)[[1]][2]) %>% unlist()

for (i in 1:length(df)){
  df[[i]]$state <- a[i]
  df[[i]]$group <- b[i]
}

df2 <- do.call('rbind',df)
df2$V7 <- if_else(df2$V7 >= 1,1,0)

df3 <-spread(data=df2, key=state,value=V7) 
df4 <-  mutate(df3,sum = rowSums(df3[,8:22]))     
c <- tapply(df4$sum, df4$group, table) 
c <- do.call('rbind',c) 
c <- data.frame(t(c))
c <- rownames_to_column(c)

d <- gather(c,key = group,value = num,2:8)
sum <- tapply(d$num, d$group, sum) %>% data.frame()
sum <- rownames_to_column(sum)
colnames(sum) <- c('group','sum')

d$sum <- sum[match(d$group,sum$group),2]
d$ratio <- d$num/d$sum
colnames(d)[1] <- 'state'
colnames(d)[2] <- 'Group'

d$Group <- if_else(d$Group == 'pqs','PQS',d$Group)
library(Hmisc)
d$Group <- capitalize(d$Group)

write.table(d,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/number_of_chromHMM_states.txt',
            sep = '\t',col.names = T,row.names = F,quote = F)

#计算state的中位数
f = c()
for (i in unique(d$Group)) {
  tmp = d[d$Group==i,]
  if (all(as.numeric(table((rep(as.numeric(tmp$state),as.numeric(tmp$num))))) == as.numeric(d[d$Group==i,"num"]))) {
    f = c(f,median(rep(as.numeric(tmp$state),as.numeric(tmp$num))))
  }
  
}
# -------------------------------------------------------------------------
#ATAC分析图
ls /home/ylxiong/G4_database/*/*/*/*/*q20.marked_duplicates.bam  > q20.marked_duplicates.bam.txt
ls /home/ylxiong/G4_database/*/*/*/*/*/*q20.marked_duplicates.bam  >> q20.marked_duplicates.bam.txt
grep ATAC q20.marked_duplicates.bam.txt > ATAC.q20.marked_duplicates.bam.txt
#合并bam文件
cd /home/ylxiong/G4_database/result/savedata2/atac_seq
samtools merge atac.bam `cat ATAC.q20.marked_duplicates.bam.txt|xargs`
qsub -q batch -V -l nodes=1:ppn=8 merge.pbs
#跨服务器复制
scp -r group/ ylxiong@211.69.141.147:/home/ylxiong/G4_database/result/savedata2/
#脚本  
cd /home/ylxiong/G4_database/result/savedata2/atac_seq
#建立索引
samtools index atac.bam atac.bam.bai
#生成 bw 文件
bamCoverage --bam atac.bam -o atac.bw --binSize 10 --normalizeUsing RPKM --extendReads 
#compute.pbs
computeMatrix reference-point -S atac.bw -R `ls ../group/*` -b 1500 -a 1500 -o matrix.mat.gz
#画图
plotProfile -m matrix.mat.gz  -out atac.pdf --colors "#FDDBC7" "#DFAFA5" "#C18383" "#A35762" "#842B40" "#67001F" "#4d4d4d"
#提交任务
qsub -q batch -V -l nodes=1:ppn=12 bw.compute.plot.pbs
# conservation-------------------------------------------------------------------------
#计算保守性分数
# phastCons 计算保守性
for i in *;do bigWigAverageOverBed /home/yulix/G4_analysis/analysis/bw_file/hg19.100way.phastCons.bw $i /home/yulix/G4_analysis/result/savedata2/conservation_score/${i/bed/phastCons};done
# phyloP  计算保守性
for i in *;do bigWigAverageOverBed /home/yulix/G4_analysis/analysis/bw_file/hg19.100way.phyloP100way.bw $i /home/yulix/G4_analysis/result/savedata2/conservation_score/${i/bed/phyloP};done

rm(list = ls());gc()
library(data.table)
library(magrittr)
library(tidyverse)

path <- "/home/yulix/G4_analysis/result/savedata2/conservation_score"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})
data <- list()
for (i in 1:length(files)){
  data[[i]] <- fread(filepath[[i]])
}
names(data) <- files
for (i in 1:length(data)){
  data[[i]]$V7 <- names(data)[i]
}

df <- do.call('rbind',data)
df$group = lapply(df$V7,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
df$V8 = lapply(df$V7,function(x)strsplit(x,'.',fixed = T)[[1]][2]) %>% unlist()

df2 <- df[,c(1,2,6,8,9)]
df3 <-spread(data=df2, key=V8,value=V6) 

colnames(df3)[1] <- 'name'
colnames(df3)[2] <- 'size'
df3 %<>% unique()
df3$group <- if_else(df3$group == 'pqs','PQS',df3$group)
library(Hmisc)
df3$group <- capitalize(df3$group)

write.table(df3,file = '/home/yulix/G4_analysis/result/savedata2/conservation_score/conservation_information',
            sep = '\t',row.names = FALSE,quote=FALSE)                  
# -------------------------------------------------------------------------
#RNA-seq
##第一步：比对
##生成hg19.splicesites.txt
python /home/ylxiong/biosoft/hisat2-2.1.0/hisat2_extract_splice_sites.py  /home/ylxiong/G4_database/data/ref/gencode.v38lift37.annotation.gtf > /home/ylxiong/G4_database/calculate/align-hisat2/hg19.splicesites.txt
qsub -q batch -V -l nodes=1:ppn=8 splice.pbs

ls /home/ylxiong/G4_database/*/*/*/*/*fq.gz > fq.txt
ls /home/ylxiong/G4_database/*/*/*/*/*/*fq.gz >> fq.txt
grep RNA fq.txt > RNA_fq.txt
#grep -v RNA fq.txt > ATAC_ChIP.fq.txt
grep Homo_sapiens RNA_fq.txt > hg_RNA_fq.txt
#grep -v Homo_sapiens RNA_fq.txt > m_RNA_fq.txt
# RNA 比对
#单端测序
for i in `less RNA_fq.single.txt2`;do
j=`basename $i`
echo "hisat2 -x /home/ylxiong/G4_database/calculate/align-hisat2/hg19/hg19 -U $i --known-splicesite-infile /home/ylxiong/G4_database/calculate/align-hisat2/hg19.splicesites.txt| samtools sort -O bam -o ${i%fq.gz}bam" > ${j%fq.gz}pbs
done
#双端测序
cat pair.txt | while read id;do arr=($id); fq1=${arr[0]}; fq2=${arr[1]};
j=`basename $fq1`
echo "hisat2 -x /home/ylxiong/G4_database/calculate/align-hisat2/hg19/hg19 -1 $fq1 -2 $fq2 --known-splicesite-infile /home/ylxiong/G4_database/calculate/align-hisat2/hg19.splicesites.txt | samtools sort -O bam  -o ${fq1%.fq.gz}.bam" > ${j%.fq.gz}.pbs
done
#提交任务
for i in *pbs;do qsub -q batch -V -l nodes=1:ppn=4 $i;done
#提取q>20的bam
ls /home/ylxiong/G4_database/*/*/*/*/*.bam > bam
ls /home/ylxiong/G4_database/*/*/*/*/*/*.bam >> bam
grep RNA bam |grep -v CHIP |grep -v Mus_musculus|grep -v q20.bam > rna.bam

cat rna.bam | while read i;do 
j=`basename $i`
echo "samtools view -b $i -q 20  -o ${i%.bam}.q20.bam" > ${j%.bam}.pbs
done
#提交任务
for i in *pbs;do qsub -q batch -V -l nodes=1:ppn=4 $i;done

#第二步：定量分析
ls /home/ylxiong/G4_database/*/*/*/*/*.q20.bam > rna.q20.bam.txt
ls /home/ylxiong/G4_database/*/*/*/*/*/*.q20.bam >> rna.q20.bam.txt
grep RNA rna.q20.bam.txt|grep -v CHIP |grep -v Mus_musculus > rna.q20.bam.txt2
grep trimmed rna.q20.bam.txt2 > rna.q20.bam.single.txt
grep -v trimmed rna.q20.bam.txt2 > rna.q20.bam.pair.txt
#单端测序
cat rna.q20.bam.single.txt | while read i;do
j=`basename $i`
echo "/home/ylxiong/biosoft/subread-2.0.3-Linux-x86_64/bin/featureCounts -t exon -g gene_id -a /home/ylxiong/G4_database/data/ref/gencode.v38lift37.annotation.gtf -o ${i%bam}counts.txt $i" > ${j%.q20.bam}.pbs
done
#双端测序
cat rna.q20.bam.pair.txt | while read i;do
j=`basename $i`
echo "/home/ylxiong/biosoft/subread-2.0.3-Linux-x86_64/bin/featureCounts -p -t exon -g gene_id -a /home/ylxiong/G4_database/data/ref/gencode.v38lift37.annotation.gtf -o ${i%bam}counts.txt $i" > ${j%.q20.bam}.pbs
done
#提交任务
for i in *pbs;do qsub -q batch -V -l nodes=1:ppn=8 $i;done

#第三步：构建表达矩阵
ls /home/ylxiong/G4_database/*/*/*/*/*.counts.txt > rna.counts.txt
ls /home/ylxiong/G4_database/*/*/*/*/*/*.counts.txt >> rna.counts.txt
grep -v Mus_musculus rna.counts.txt > hg.rna.counts.txt
grep Mus_musculus rna.counts.txt > m.rna.counts.txt

for i in `less hg.rna.counts.txt`;do j=`basename $i`
sed '1,2d' $i|cut -f7 > ${j%.q20.counts.txt}.a;done
ls *.a |xargs -n 1 > gene_id
paste `ls *.a |xargs` -d "\t" > counts.matrix
sed '1,2d' SRR11113121_trimmed.q20.counts.txt |awk '{print $1}' > Geneid
paste Geneid counts.matrix -d "\t" > counts.matrix2

#计算 tpm
rm(list = ls());gc();rm(list = ls())
counts <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/counts.matrix2') %>% data.frame()
a <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/gene_id',header = F) 
a$V1 %<>% lapply(.,function(x)strsplit(x,'_',fixed = T)[[1]][1]) %>% unlist()
colnames(counts)[2:42] <- a$V1

exon_length <- fread('/home/yulix/G4_analysis/data/ref/Human.nonreExon_length.bed') %>% data.frame()
counts$V1 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
gene.length <- tapply(exon_length$V5, exon_length$V4, sum) %>% as.data.frame()
# tpm
tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  row.names(tpm) = row.names(exprset)
  return(tpm)
}
counts$length = gene.length[match(counts$V1,row.names(gene.length)),1]
tpm = tpm.calculate(counts[,-c(1,ncol(counts))],counts$length)

gene_exp <- data.frame(counts$V1,tpm)
colnames(gene_exp)[1] <- 'gene_id'

gene_list <- read.csv('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/gene_list3',
                      sep = '\t',header = F,stringsAsFactors = F)
gene_list$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
gene_exp$chr <- gene_list[match(gene_exp$gene_id,gene_list$V4),1]
gene_exp$start <- gene_list[match(gene_exp$gene_id,gene_list$V4),2]
gene_exp$end <- gene_list[match(gene_exp$gene_id,gene_list$V4),3]
gene_exp <- gene_exp[,c(43:45,1:42)]
write.table(gene_exp,file = '/home/yulix/G4_analysis/result/savedata2/RNA_seq/gene_exp.bed',
            sep = "\t",quote = F,row.names = F,col.names = T)


#group2.bed 的数据有问题，进行修改
awk '{if($7=="Group2") print $0}'  6groups.bed |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' |bedtools sort -i > group2.bed

rm(list = ls());gc();rm(list = ls())
G4 <- fread('~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed') %>% data.frame()
groups <- fread('~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed') %>% data.frame()
G4$V7 <-'G4'
data <- rbind(G4,groups)
write.table(data,file = '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4_group.bed',
            sep = "\t",quote = F,row.names = F,col.names = F)
#fraction.G4_group------------------------------------------------------------------------
for i in /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*;do j=`basename $i`
bedtools intersect -a SRR11113121.bed -b $i -c > ./SRR11113121_c/${j%bed}txt;done

path <- '/home/yulix/G4_analysis/result/savedata2/RNA_seq/SRR11113121_c'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

a <- c(paste0('Group',1:6))
for (i in 1:length(data)){
  data[[i]]$group <- a[i]
}

data2 <-do.call('rbind',data)
data3 <- data2[data2$V6 != 0,]
data4 <- tapply((data3$V5 > 1), data3$group, table) 

for (i in 1:length(data)){
  data4[[i]] %<>% data.frame()
  data4[[i]]$group <- names(data4)[[i]]
}

df5 <-do.call('rbind',data4)
df6<-spread(data=df5,key=Var1,value=Freq)
df6$ratio <- df6$`TRUE`/(df6$`TRUE`+df6$`FALSE`)     
df6 <- df6[,c(1,4)]      
colnames(df6) <- c('Class','Fraction')
  
G4.num <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/G4.num.txt') %>% data.frame()
table(G4.num$class)
#G4-containing genes   G4-depleted genes 
#25322               37124 
nrow(G4.num[G4.num$class == 'G4-depleted genes' & G4.num$expression > 1,])/table(G4.num$class)[[2]]  #0.1974195
nrow(G4.num[G4.num$class == 'G4-containing genes' & G4.num$expression > 1,])/table(G4.num$class)[[1]] #0.532067

df <- data.frame(a = 0.1974195,
                 b = 0.532067)
colnames(df) <- c('G4-depleted genes','G4-containing genes')
df <-  t(as.matrix(df)) %>% data.frame()
df <- rownames_to_column(df)
colnames(df) <- c('Class','Fraction')

df <- rbind(df,df6)
write.table(df,file = '/home/yulix/G4_analysis/result/savedata2/RNA_seq/fraction.G4_group.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)
# -------------------------------------------------------------------------
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' gene_exp.bed |sed '1d' > SRR11113121.bed
bedtools intersect -a SRR11113121.bed -b /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed -c > SRR11113121_G4
bedtools intersect -a SRR11113121.bed -b /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed -wo > SRR11113121_G4_length

# expression_level.G4_group-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
G4.num <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/G4.num.txt') %>% data.frame()
G4.num$class %<>% factor(.,levels = c('G4-depleted genes','G4-containing genes'))

path <- '/home/yulix/G4_analysis/result/savedata2/RNA_seq/SRR11113121_c'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

a <- c(paste0('Group',1:6))
for (i in 1:length(data)){
  data[[i]]$class <- a[i]
}
data2 <-do.call('rbind',data)
data3 <- data2[data2$V6 != 0,]
colnames(data3)[5] <- 'expression'

df <- rbind(data3,G4.num)
write.table(df,file = '/home/yulix/G4_analysis/result/savedata2/RNA_seq/expression_level.G4_group.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# promoter_fraction.G4_group-------------------------------------------------------------------------
for i in /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*;do j=`basename $i`
bedtools intersect -a /home/yulix/G4_analysis/result/savedata/epigenetics/genome_region/promoter2.bed -b $i -c > ./promoter_SRR11113121_c/${j%bed}txt;done

bedtools intersect -a /home/yulix/G4_analysis/result/savedata/epigenetics/genome_region/promoter2.bed -b /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed -c > ./promoter_SRR11113121_c/G4.txt

rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata2/RNA_seq/promoter_SRR11113121_c'
files <- dir(path)[-1]
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

a <- c(paste0('Group',1:6))
for (i in 1:length(data)){
  data[[i]]$group <- a[i]
}

data2 <-do.call('rbind',data)
data2 <- data2[data2$V6 != 0,]

g4 <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/promoter_SRR11113121_c/G4.txt') %>% data.frame()
g4$group <- if_else(g4$V6 ==0,'G4-depleted genes','G4-containing genes')

df <- rbind(g4,data2)

df$V4 %<>%  lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
exp <- fread('~/G4_analysis/result/savedata2/RNA_seq/SRR11113121.bed') %>% data.frame()
df$exp <- exp[match(df$V4,exp$V4),5]

write.table(df,file = '/home/yulix/G4_analysis/result/savedata2/RNA_seq/promoter_expression_level.G4_group.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

data4 <- tapply((df$exp > 1), df$group, table) 

for (i in 1:length(data4)){
  data4[[i]] %<>% data.frame()
  data4[[i]]$group <- names(data4)[[i]]
}

df5 <-do.call('rbind',data4)
df6<-spread(data=df5,key=Var1,value=Freq)
df6$ratio <- df6$`TRUE`/(df6$`TRUE`+df6$`FALSE`)     
df6 <- df6[,c(1,4)]      
colnames(df6) <- c('Class','Fraction')

write.table(df6,file = '/home/yulix/G4_analysis/result/savedata2/RNA_seq/promoter_fraction.G4_group.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# enhancer-------------------------------------------------------------------------
# pqs,G4 和6个group promoter 区以外的部分
for i in *;do bedtools subtract -a $i -b /home/ylxiong/G4_database/result/savedata2/epigenetics/promoter2.bed > /home/ylxiong/G4_database/result/savedata2/epigenetics/group_subtract_e/${i%bed}promoter.bed;done

for i in *;do bedtools intersect -a $i -b /home/sqian/hs.bed -c > /home/ylxiong/G4_database/result/savedata2/epigenetics/enhancer/enhancer_${i%.bed};done

awk '{print $3"\t"$4"\t"$5}' /home/sqian/SE_package.bed |sed '1d' > SE_package.bed
for i in *;do bedtools intersect -a $i -b /home/ylxiong/G4_database/result/savedata2/epigenetics/SE_package.bed -c > /home/ylxiong/G4_database/result/savedata2/epigenetics/enhancer/Senhancer_${i%.bed};done
#结果复制到实验室服务器
scp -r ylxiong@211.69.141.147:/home/ylxiong/G4_database/result/savedata2/epigenetics/enhancer/ ./

rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata/epigenetics/enhancer'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

res <- vector()
for (i in 1:length(data)){
  res[i] <- nrow(data[[i]][data[[i]]$V7 != 0,])/nrow(data[[i]])
}

res <- data.frame(ratio = res,V2 = files)
res <- separate(res,col = V2,sep = '_',into = c('class','group'))
res$group %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
res$group <- if_else(res$group == 'pqs','PQS',res$group)
library(Hmisc)
res$group <- capitalize(res$group)

write.table(res,file = '/home/yulix/G4_analysis/result/savedata/epigenetics/enhancer.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# -------------------------------------------------------------------------
#处理成规范的bed文件格式
#for i in 3_UTR.bed 5_UTR.bed exon.bed ;do awk '{print $1"\t"$2"\t"$3"\t"NR"\t""1""\t"$4}' $i > ${i%.bed}1.bed;done

# 相同链上 4个基因区域分别和 PQS ,G4 ,6 groups的相交数量
for i in 3_UTR1.bed 5_UTR1.bed exon1.bed intron.bed;do
for j in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;
do k=`basename $j`
bedtools intersect -a $i -b $j -c -s > ./num/${i%.bed}_${k%.bed}.s
done
done
# 相反链上 4个基因区域分别和 PQS ,G4 ,6 groups 的相交数量
for i in 3_UTR1.bed 5_UTR1.bed exon1.bed intron.bed;do
for j in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;
do k=`basename $j`
bedtools intersect -a $i -b $j -c -S > ./num/${i%.bed}_${k%.bed}.S
done
done

# 相同链上 4个基因区域分别和 PQS ,G4 ,6 groups的相交长度
for i in 3_UTR1.bed 5_UTR1.bed exon1.bed intron.bed;do
for j in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;
do k=`basename $j`
bedtools intersect -a $i -b $j -wo -s > ./length/${i%.bed}_${k%.bed}.s
done
done
# 相反链上 4个基因区域分别和 PQS ,G4 ,6 groups 的相交长度
for i in 3_UTR1.bed 5_UTR1.bed exon1.bed intron.bed;do
for j in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;
do k=`basename $j`
bedtools intersect -a $i -b $j -wo -S > ./length/${i%.bed}_${k%.bed}.S
done
done

# overlap_length.strand-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata/gene_region_intersect_G4/length'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)data.table::fread(x))

names(data) <- files
for (i in 1:length(data)){
  data[[i]]$V14 = names(data)[i]
}

df = do.call('rbind',data)
df$length <- df$V3-df$V2
df <- df[,13:15]
df <- df[df$length !=0,]
df$density <- df$V13/df$length


df$V5 <- lapply(df$V14,function(x)strsplit(x,'.',fixed = T)[[1]][2]) %>% unlist()
df <- mutate(df,strand = if_else(df$V5 == 's','same_strandedness','opposite_strandedness'))

df$class <- lapply(df$V14,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
df$class %<>% gsub('3_UTR1','3\'UTR',.)  %>% gsub('5_UTR1','5\'UTR',.) %>% gsub('exon1','Exon',.)

df <- separate(df,class,sep = '_',into = c('class','type'))              
df$class <- Hmisc::capitalize(df$class)

df$type %<>% gsub('pqs','PQS',.)
df$type <- Hmisc::capitalize(df$type)

df2 <- df[,c(4,6:8)]

df3 <- aggregate(df2$density, by = list(df2$class,df2$type,df2$strand), mean) %>% data.frame()

colnames(df3) <- c('class','type','strand','mean')
write.table(df3,file = '/home/yulix/G4_analysis/result/savedata/gene_region_intersect_G4/overlap_length.strand.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# overlap_num.strand-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata/gene_region_intersect_G4/num'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)data.table::fread(x))

names(data) <- files
for (i in 1:length(data)){
  data[[i]]$V14 = names(data)[i]
}

df = do.call('rbind',data)
df$length <- df$V3-df$V2
df <- df[,7:9]
df <- df[df$length !=0,]
df$density <- df$V7/df$length


df$V5 <- lapply(df$V14,function(x)strsplit(x,'.',fixed = T)[[1]][2]) %>% unlist()
df <- mutate(df,strand = if_else(df$V5 == 's','same_strandedness','opposite_strandedness'))

df$class <- lapply(df$V14,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
df$class %<>% gsub('3_UTR1','3\'UTR',.)  %>% gsub('5_UTR1','5\'UTR',.) %>% gsub('exon1','Exon',.)

df <- separate(df,class,sep = '_',into = c('class','type'))              
df$class <- Hmisc::capitalize(df$class)

df$type %<>% gsub('pqs','non-eG4',.)
df$type <- Hmisc::capitalize(df$type)

df2 <- df[,c(4,6:8)]

df3 <- aggregate(df2$density, by = list(df2$class,df2$type,df2$strand), mean) %>% data.frame()

colnames(df3) <- c('class','type','strand','mean')
df3$mean <- df3$mean*(10^4)

write.table(df3,file = '/home/yulix/G4_analysis/result/savedata/gene_region_intersect_G4/overlap_num.strand.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)
# -------------------------------------------------------------------------
rm(list = ls());gc()
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/repeats_hg19_ucsc.sort') %>% data.table()
a = c('DNA','LINE','LTR','SINE')
data <- list()
for (i in 1:length(a)){
  data[[i]] <- df[df$V7 == a[i],]
}

df <- do.call('rbind',data)
write.table(df,file = '/home/yulix/G4_analysis/result/savedata/G4_TE/DNA_LINE_LTR_SINE.txt',
            sep = "\t",quote = F,row.names = F,col.names = F)

#合成一个总表
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""PQS"}' pqs.bed > pqs.2.bed
cat pqs.2.bed G4_group.bed |bedtools sort -i - > G4_group_PQS.bed

#对PQS随机排序 10次
cd ~/G4_analysis/result/savedata/G4_TE/random_PQS
for i in {1..10};do sort --random-sort ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.2.bed |head -391404 | bedtools sort -i - > PQS.${i};done

cd ~/G4_analysis/result/savedata/G4_TE
for i in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;do
j=`basename $i`
bedtools intersect -a DNA_LINE_LTR_SINE.txt -b $i -c > ./G4_intersect_TE/${j%.bed}
done

cd ~/G4_analysis/result/savedata/G4_TE/random_PQS
for i in *;do
bedtools intersect -a ../DNA_LINE_LTR_SINE.txt -b $i -c > ../G4_intersect_TE/${i}
done

#shuffle 10 次
for i in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;do
k=`basename $i`
for j in {1..10}
do bedtools shuffle -i $i -g my.genome -noOverlapping -chrom -seed $j > ./shuffle_files/${k%.bed}_${j}
done
done

# PQS shuffle 10 次
sort --random-sort pqs.2.bed |head -391404 | bedtools sort -i - > pqs.some.bed
cd ~/G4_analysis/result/savedata/G4_TE/shuffle_files
for i in {1..10};
do bedtools shuffle -i /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.some.bed -g ../my.genome -noOverlapping -chrom -seed $i > PQS.${i};done

cd ~/G4_analysis/result/savedata/G4_TE/shuffle_files
for i in *;do
bedtools intersect -a ../DNA_LINE_LTR_SINE.txt -b $i -c > ../shuffle_G4_intersect_TE/${i}
done

# fold_enrichment.PQS_G4-------------------------------------------------------------------------

rm(list = ls());gc();rm(list = ls())
library(data.table)
library(magrittr)
library(tidyverse)
path <- '/home/yulix/G4_analysis/result/savedata/G4_TE/G4_intersect_TE'
files <- dir(path)
filespath <- sapply(files, function(x)paste(path,x,sep = '/'))

obs_g4 <- list()
obs_g4 <- lapply(filespath, function(x)fread(x))

for (i in 1:length(obs_g4)){
  obs_g4[[i]] = obs_g4[[i]][obs_g4[[i]]$V8 != 0,]
}

a <- vector()
for (i in 1:length(obs_g4)){
  a[i] = nrow(obs_g4[[i]])
}

df <- data.frame(obs = a,group = names(obs_g4))
df$group %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
df.1 <- aggregate(df$obs,by = list(df$group),mean)
colnames(df.1) <- c('group','obs')

path2 <- '/home/yulix/G4_analysis/result/savedata/G4_TE/shuffle_G4_intersect_TE'
files2 <- dir(path2)
filespath2 <- sapply(files2, function(x)paste(path2,x,sep = '/'))

simu <- list()
simu <- lapply(filespath2, function(x)fread(x)) 

for (i in 1:length(simu)){
  simu[[i]] = simu[[i]][simu[[i]]$V8 != 0,]
}

b <- vector()
for (i in 1:length(simu)){
  b[i] = nrow(simu[[i]])
}

df2 <- data.frame(simu = b,group = names(simu))
df2$group %<>% lapply(.,function(x)strsplit(x,'_',fixed = T)[[1]][1]) %>% unlist
df2$group %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist
df2.1 <- aggregate(df2$simu,by = list(df2$group),mean)
colnames(df2.1) <- c('group','simu')

res <- left_join(df.1,df2.1,by = c('group' = 'group'))
res$FoldEnrichment <- res$obs/res$simu
res$group <- Hmisc::capitalize(res$group)

write.table(res,file = '/home/yulix/G4_analysis/result/savedata/G4_TE/fold_enrichment.PQS_G4.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# -------------------------------------------------------------------------
for i in `ls /home/yulix/G4_analysis/result/savedata/G4_TE/random_PQS/*` `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/group*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;do j=`basename $i`
bedtools intersect -a gene_list3 -b $i -c > ./intersect_obs/${j};done

for i in `ls ~/G4_analysis/result/savedata/G4_TE/shuffle_files/*`;do j=`basename $i`
bedtools intersect -a gene_list3 -b $i -c > ./intersect_simu/${j};done

rm(list = ls());gc();rm(list = ls())
phenotypic <- read.csv('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_merge',
                       sep = '\t',header = F,stringsAsFactors = F)
tsg <- read.csv('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/tsg_merge.1',
                sep = '\t',header = F,stringsAsFactors = F) 
oncogene <- read.csv('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/oncogene_merge.1',
                     sep = '\t',header = F,stringsAsFactors = F)

phenotypic <- phenotypic[,c(1,3)]
phenotypic$V3 %<>% gsub("[0-9]","",.)
colnames(phenotypic)[2] <- 'class'
tsg$class <- 'tsg'
oncogene$class <- 'oncogene'
df <- rbind(phenotypic,tsg,oncogene)
colnames(df)[1] <- 'gene_id' 
df$class <- Hmisc::capitalize(df$class)

write.table(df,file = '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_tsg_oncogene.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# fold_enrichment.desease-------------------------------------------------------------------------

rm(list = ls());gc();rm(list = ls())
library(data.table)
library(magrittr)
library(tidyverse)
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_tsg_oncogene.txt') %>% data.frame()

path <- '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/intersect_obs'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))


obs <-data[[1]][,1:4] %>% as.data.frame()

for (i in 1:length(data)){
  obs[,files[i]] = data[[i]][,9]
}
obs$PQS.mean <- apply(obs[,12:21],1,mean) 
obs <- obs[,c(1:11,22)]
obs$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist() 

df[,3:10] = obs[match(df$gene_id,obs$V4),5:12] 

sum(is.na(df))
df <- na.omit(df)
colnames(df)[3:10] <- c('G4',paste0('Group',1:6),'PQS') 

sum_group <- lapply(df[,3:10], function(x)tapply(x,df$class,sum))
sum_group.1 = lapply(sum_group,function(x)x%>% unlist() %>% data.frame() %>% rownames_to_column(.))

for (i in 1:length(sum_group.1)){
  sum_group.1[[i]]$group <- names(sum_group)[i]
}

sum_group_obs <- do.call('rbind',sum_group.1)
colnames(sum_group_obs) <- c('class','obs','group')

path.1 <- '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/intersect_simu'
files.1 <- dir(path.1)
filespath.1 <- lapply(files.1, function(x)paste(path.1,x,sep = '/'))
data.1 <- list()
data.1 <- lapply(filespath.1,function(x)fread(x))

simu <-data.1[[1]][,1:4] %>% as.data.frame()

for (i in 1:length(data.1)){
  simu[,files.1[i]] = data.1[[i]][,9]
}

simu[85:92] <- lapply(seq(5,75,10),function(x)apply(simu[,x:(x+9)],1,mean))
simu.1 <- simu[,c(1:4,85:92)]
simu.1$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist() 
colnames(simu.1)[5:12] <- c('G4',paste0('Group',1:6),'PQS') 

df.1 = fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_tsg_oncogene.txt') %>% data.frame()
df.1[,3:10] = simu.1[match(df.1$gene_id,simu.1$V4),5:12] 

sum(is.na(df.1))
df.1 <- na.omit(df.1)

sum_group2 <- lapply(df.1[,3:10], function(x)tapply(x,df.1$class,sum))
sum_group2.1 = lapply(sum_group2,function(x)x%>% unlist() %>% data.frame() %>% rownames_to_column(.))

for (i in 1:length(sum_group2.1)){
  sum_group2.1[[i]]$group <- names(sum_group2.1)[i]
}

sum_group_simu <- do.call('rbind',sum_group2.1)
colnames(sum_group_simu) <- c('class','simu','group')

#a = cbind(sum_group_obs,sum_group_simu$simu)
sum_group_obs$`Fold enrichment` <- sum_group_obs$obs/sum_group_simu$simu
sum_group_obs$class %<>% gsub('Tsg','TSG',.)

write.table(sum_group_obs,file = '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/fold_enrichment.desease.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)

# promoter_fold_enrichment.desease-------------------------------------------------------------------------
for i in `ls /home/yulix/G4_analysis/result/savedata/G4_TE/random_PQS/*` `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/group*` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;do j=`basename $i`
bedtools intersect -a /home/yulix/G4_analysis/result/savedata/epigenetics/genome_region/promoter2.bed -b $i -c > ./promoter_obs/${j};done

for i in `ls ~/G4_analysis/result/savedata/G4_TE/shuffle_files/*`;do j=`basename $i`
bedtools intersect -a /home/yulix/G4_analysis/result/savedata/epigenetics/genome_region/promoter2.bed -b $i -c > ./promoter_simu/${j};done

rm(list = ls());gc();rm(list = ls())
library(data.table)
library(magrittr)
library(tidyverse)
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_tsg_oncogene.txt') %>% data.frame()

path <- '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/promoter_obs'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))


obs <-data[[1]][,1:4] %>% as.data.frame()

for (i in 1:length(data)){
  obs[,files[i]] = data[[i]][,6]
}
obs$PQS.mean <- apply(obs[,12:21],1,mean) 
obs <- obs[,c(1:11,22)]
obs$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist() 

df[,3:10] = obs[match(df$gene_id,obs$V4),5:12] 

sum(is.na(df))
df <- na.omit(df)
colnames(df)[3:10] <- c('G4',paste0('Group',1:6),'PQS') 

sum_group <- lapply(df[,3:10], function(x)tapply(x,df$class,sum))
sum_group.1 = lapply(sum_group,function(x)x%>% unlist() %>% data.frame() %>% rownames_to_column(.))

for (i in 1:length(sum_group.1)){
  sum_group.1[[i]]$group <- names(sum_group)[i]
}

sum_group_obs <- do.call('rbind',sum_group.1)
colnames(sum_group_obs) <- c('class','obs','group')

path.1 <- '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/promoter_simu'
files.1 <- dir(path.1)
filespath.1 <- lapply(files.1, function(x)paste(path.1,x,sep = '/'))
data.1 <- list()
data.1 <- lapply(filespath.1,function(x)fread(x))

simu <-data.1[[1]][,1:4] %>% as.data.frame()

for (i in 1:length(data.1)){
  simu[,files.1[i]] = data.1[[i]][,6]
}

simu[85:92] <- lapply(seq(5,75,10),function(x)apply(simu[,x:(x+9)],1,mean))
simu.1 <- simu[,c(1:4,85:92)]
simu.1$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist() 
colnames(simu.1)[5:12] <- c('G4',paste0('Group',1:6),'PQS') 

df.1 = fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/phenotypic_tsg_oncogene.txt') %>% data.frame()
df.1[,3:10] = simu.1[match(df.1$gene_id,simu.1$V4),5:12] 

sum(is.na(df.1))
df.1 <- na.omit(df.1)

sum_group2 <- lapply(df.1[,3:10], function(x)tapply(x,df.1$class,sum))
sum_group2.1 = lapply(sum_group2,function(x)x%>% unlist() %>% data.frame() %>% rownames_to_column(.))

for (i in 1:length(sum_group2.1)){
  sum_group2.1[[i]]$group <- names(sum_group2.1)[i]
}

sum_group_simu <- do.call('rbind',sum_group2.1)
colnames(sum_group_simu) <- c('class','simu','group')

sum_group_obs$`Fold enrichment` <- sum_group_obs$obs/sum_group_simu$simu
sum_group_obs$class %<>% gsub('Tsg','TSG',.)

write.table(sum_group_obs,file = '/home/yulix/G4_analysis/result/savedata/G4_TE/disease/promoter_fold_enrichment.desease.txt',
            sep = "\t",quote = F,row.names = F,col.names = T)
# -------------------------------------------------------------------------
#鸡基因组run pqsfinder
samtools faidx Gallus_gallus.genome.fa

rm(list = ls());gc();rm(list = ls())
library(pqsfinder)
library(Biostrings)
genome <- readDNAStringSet("/home/yulix/G4_analysis/data/ref/Gallus_gallus/Gallus_gallus.genome.fa") 

chr_pqs <- list()
for (i in 1:34){
  chr_pqs[[i]] <- pqsfinder(genome[[i]],overlapping = FALSE, min_score = 50)
}

a <- c(10,11,12,13,14,15,16,17,18,19,1,20,21,22,23,24,25,26,27,28,2,30,31,32,33,3,4,5,6,7,8,9,'W','Z')

a <- c(paste0('chr',a))
list <- list()
for (i in 1:length(chr_pqs)){
  list[[i]] <- data.frame(list(chr_pqs[[i]]@ranges,chr_pqs[[i]]@elementMetadata))
  list[[i]]$chr <- a[i]
}

df <- do.call('rbind',list)
df <- df[,c(15,1:14)]
write.table(df,file="/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/gallus.pqs.txt",
            sep = '\t',row.names = F,col.names = F,quote=FALSE)

# -------------------------------------------------------------------------
awk '{print $1"\t"$2"\t"$3"\t""m_"NR"\t"$6"\t"$5}' m.pqs.txt > m.pqs.txt2
for i in *narrowPeak;do bedtools intersect -a /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/m.pqs.txt2 -b $i -c > ${i%narrowPeak}c;done

# count.matrix 添加列名-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
df <- fread('/home/yulix/G4_analysis/result/savedata2/ligand/counts.matrix2') %>% data.frame()
df2 <- read.table('/home/yulix/G4_analysis/result/savedata2/ligand/sample.txt', 
                  sep = '\t',header = F,stringsAsFactors = F)
colnames(df)[1] <- 'gene_id'
colnames(df)[2:ncol(df)] <- df2$V2

write.table(df,file = '/home/yulix/G4_analysis/result/savedata2/ligand/counts.colnames.matrix',
            sep = "\t",quote = F,row.names = F,col.names = T)

# 差异表达分析deseq2-------------------------------------------------------------------------
#表达矩阵
#剔除Y上的同源基因
grep -v Y counts.colnames.matrix > counts.colnames.matrix2

rm(list = ls());gc();rm(list = ls())
library(magrittr)
gene_exp <- read.csv("/home/yulix/G4_analysis/result/savedata2/ligand/counts.colnames.matrix2",sep = "")
gene_exp$gene_id = gsub("\\..*","",gene_exp$gene_id)
rownames(gene_exp) <- gene_exp$gene_id
gene_exp <- gene_exp[,-1]
##total
dim(gene_exp) # 62397    18
#样本信息表
colnames(gene_exp)
sample_info <- data.frame(group = colnames(gene_exp))
rownames(sample_info) <- colnames(gene_exp)
sample_info$type <- lapply(sample_info$group,function(x)strsplit(x,'_',fixed = T)[[1]][1]) %>% unlist()
##差异分析
library(DESeq2)

sample_info$type %<>% factor()
sample_info$type <- relevel(sample_info$type, ref = "Untreated")
dds <- DESeqDataSetFromMatrix(countData = gene_exp, colData = sample_info, design = ~ type)
dds <- DESeq(dds) 
resultsNames(dds)
res <- list()
for (i in resultsNames(dds)[-1]){
  res[i] <- results (dds,name = i)
}

##差异表达基因
res2 <- lapply(res, function(x)x[!is.na(x$padj) & x$padj<=0.05,])
##上调基因
res3<-lapply(res2,function(x)x[x$log2FoldChange>=1,])
##下调基因
res4<-lapply(res2,function(x)x[x$log2FoldChange<=-1,])

up <- list()
down <- list()
for (i in 3:5){
  up[[i]] = data.frame(gene_id = res3[[i]]@rownames,log2FoldChange = res3[[i]]@listData$log2FoldChange)
  down[[i]] = data.frame(gene_id = res4[[i]]@rownames,log2FoldChange = res4[[i]]@listData$log2FoldChange)
  up[[i]]$class <- names(res3)[i] 
  down[[i]]$class <- names(res4)[i] 
}

up.df <- do.call('rbind',up[3:5])
up.df$direction <- 'up'
down.df <- do.call('rbind',down[3:5])  
down.df$direction <- 'down'
resdata <- rbind(up.df,down.df)
# 得到csv格式的差异表达分析结果
write.csv(resdata,file= "/home/yulix/G4_analysis/result/savedata2/ligand/up_down_gene.txt",row.names = F)

gene_list <- read.csv('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/gene_list3',
                      sep = '\t',header = F,stringsAsFactors = F)
gene_list$V4 %<>% lapply(.,function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()
resdata[5:7] <- gene_list[match(resdata$gene_id,gene_list$V4),1:3]
resdata <- resdata[,c(5:7,1:4)]
write.table(resdata,file = '/home/yulix/G4_analysis/result/savedata2/ligand/up_down_gene.bed',
            sep = '\t',col.names = F,row.names = F,quote = F)
# -------------------------------------------------------------------------
for i in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*`;do j=`basename $i`
bedtools intersect -a up_down_gene.bed -b $i -c > ./up_down_gene.G4/${j%.bed};done

rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata2/ligand/up_down_gene.G4/'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)data.table::fread(x))
names(data) <- files

for (i in 1:length(data)){
  data[[i]]$group <- names(data)[i]
}
df <- do.call('rbind',data)

df.1 <- aggregate(df$V8,by = list(df$V6,df$V7,df$group),sum)
df.2 <- data.frame(num = c(123053,75922,52136,59175,53742,27376,1118855),
                   group = c(paste0('Group',1:6),'PQS'))

df.1$Group.3 %<>% gsub('pqs','PQS',.) %>% gsub('group','Group',.)
df.1$sum <- df.2[match(df.1$Group.3,df.2$group),1]
df.1$ratio <- df.1$x/df.1$sum
df.1$Group.4 <- lapply(df.1$Group.1,function(x)strsplit(x,'_',fixed = T)[[1]][2]) %>% unlist()

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
library(ggplot2)
ggplot(data = df.1,aes(x = Group.4,y= ratio,fill = Group.3)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab('Fraction of overlapped G4 in each group') +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,0.08))+
  facet_wrap(.~Group.2)
# -------------------------------------------------------------------------
awk '$6=="type_PDS_vs_Untreated" || $7=="up" {print $0}' up_down_gene.bed > up_down_gene.class/PDS.up.txt
awk '$6=="type_PhenDC3_vs_Untreated" || $7=="up" {print $0}' up_down_gene.bed > up_down_gene.class/PhenDC3.up.txt
awk '$6=="type_PhenDC6_vs_Untreated" || $7=="up" {print $0}' up_down_gene.bed > up_down_gene.class/PhenDC6.up.txt
awk '$6=="type_PDS_vs_Untreated" || $7=="down" {print $0}' up_down_gene.bed > up_down_gene.class/PDS.down.txt
awk '$6=="type_PhenDC3_vs_Untreated" || $7=="down" {print $0}' up_down_gene.bed > up_down_gene.class/PhenDC3.down.txt
awk '$6=="type_PhenDC6_vs_Untreated" || $7=="down" {print $0}' up_down_gene.bed > up_down_gene.class/PhenDC6.down.txt

for a in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*`;do 
for b in *;do j=`basename $a`;k=`basename $b`
bedtools intersect -a $a -b $b -c > ../up_down_gene.G4/${j%.bed}.${k%.txt};done;done

rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata2/ligand/up_down_gene.G4'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)data.table::fread(x))
names(data) <- files

for (i in 1:length(data)){
  data[[i]]$class <- names(data)[i]
}
df <- do.call('rbind',data)
df %<>% separate(class, c('group','condition','direction'))

df.1 <- aggregate(df$V7,by = list(df$group,df$condition,df$direction),sum)
df.2 <- data.frame(num = c(3747,4342,3711,4266,2718,3488),
                   group = c("PDS.down.txt","PDS.up.txt","PhenDC3.down.txt","PhenDC3.up.txt","PhenDC6.down.txt","PhenDC6.up.txt"))
df.2$group %<>% gsub('.txt','',.)
df.2 %<>% separate(group, c('condition','direction'))

df.1$sum <- df.2[match(paste0(df.1$Group.2,df.1$Group.3),paste0(df.2$condition,df.2$direction)),1]
df.1$ratio <- df.1$x/df.1$sum
df.1$Group.1 %<>% gsub('pqs','PQS',.) %>% gsub('g','G',.)

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
library(ggplot2)
ggplot(data = df.1,aes(x = Group.2,y= ratio,fill = Group.1)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab('Fraction of overlapped gene in each group') +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,25))+
  facet_wrap(.~Group.3)
# -------------------------------------------------------------------------
bedtools intersect -a SRR11113121.bed -b ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed -c > ./SRR11113121_c/pqs.txt

rm(list = ls());gc();rm(list = ls())
df <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/expression_level.G4_group.txt') %>% data.frame()
df2 <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/SRR11113121_c/pqs.txt') %>% data.frame()
df2 <- df2[df2$V6 != 0,]
df2$class <- 'PQS'
colnames(df)[5] <- 'V5'
df3 <- rbind(df,df2)
df3 <- df3[,c(1:4,7)]
df3 <- unique(df3)
write.table(df3,file = '/home/yulix/G4_analysis/result/savedata2/ligand/G4_contain_gene.bed',
            sep = '\t',col.names = F,row.names = F,quote = F)

for i in *;do j=`basename $i`
bedtools intersect -a /home/yulix/G4_analysis/result/savedata2/ligand/G4_contain_gene.bed -b $i -c > ../up_down_gene.G4/${j%.txt};done

rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata2/ligand/up_down_gene.G4'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)data.table::fread(x))
names(data) <- files

for (i in 1:length(data)){
  data[[i]]$class <- names(data)[i]
}
df <- do.call('rbind',data)
df %<>% separate(class, c('condition','direction'))

df$V9 <- 1
df.1 <- aggregate(df$V9,by = list(df$V5,df$condition,df$direction),sum)

df <- df[df$V6 != 0,] 
df.2 <- aggregate(df$V9,by = list(df$V5,df$condition,df$direction),sum)

df.1$ratio <- df.2$x/df.1$x

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('black','#4d4d4d','#b2182b',b)
library(ggplot2)
df.1$Group.1 %<>% factor(.,levels = c('G4-depleted genes','PQS','G4-containing genes',
                                      paste0('Group',1:6)))

ggplot(data = df.1,aes(x = Group.2,y= ratio,fill = Group.1)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab('Fraction of overlapped gene in each group') +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "top" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,0.2))+
  facet_wrap(.~Group.3)

# GO-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
path <- "/home/yulix/G4_analysis/result/savedata2/ligand/up_down_gene.class"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

a = files %>% gsub('.txt','',.)

peak=list(PDS.down= filepath[[1]],PDS.up = filepath[[2]],PhenDC3.down = filepath[[3]],
          PhenDC3.up = filepath[[4]],PhenDC6.down = filepath[[5]],PhenDC6.up = filepath[[6]])
          
peakAnnoList = list()
for (i in 1:length(peak)) {
  peakAnnoList[[i]] <- lapply(peak[[i]], annotatePeak, TxDb=txdb,
                              tssRegion=c(-1500, 1500), verbose=FALSE)
}        
gene=list()
gene= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
gene = lapply(gene, function(x)data.frame(x))

down <- intersect(gene[[1]]$x,gene[[3]]$x)
down <- intersect(down,gene[[5]]$x)

up <- intersect(gene[[2]]$x,gene[[4]]$x)
up <- intersect(up,gene[[6]]$x)

erich.go.BP = enrichGO(gene = down,
                              OrgDb = org.Hs.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1)

erich.go.BP@result$Description %<>% Hmisc::capitalize(.)
ggplot(erich.go.BP, aes(GeneRatio, Description), showCategory=10) +
  geom_point(aes(color=p.adjust, size=Count)) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text= element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_blank())

erich.go.BP.up = enrichGO(gene = up,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1)

erich.go.BP.up@result$Description %<>% Hmisc::capitalize(.)
ggplot(erich.go.BP.up, aes(GeneRatio,Description), showCategory=10) +
  geom_point(aes(color=p.adjust, size=Count)) +
  scale_color_continuous(low='#fddbc7', high='#b2182b') +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text= element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_blank())
# TF-binding diversity-------------------------------------------------------------------------
for a in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*`;do j=`basename $a`
bedtools intersect -a $a -b Human.macs2.clusters.interval.merge3.hg19.bed2 -wa -wb > ./G4.TF/${j%.bed};done

for i in *;do awk 'BEGIN{OFS="\t"} {print $4,$10}' $i |sort |uniq > $i.uniq;done

rm(list = ls());gc();rm(list = ls())
path <- "~/G4_analysis/result/savedata2/TF/G4.TF"
files <- dir(path)[seq(2,14,2)]
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

a = files %>% gsub('.uniq','',.) %>% gsub('pqs','PQS',.)                   
a = Hmisc::capitalize(a)

data <- list()
data <- lapply(filepath,function(x)read.csv(x,sep = '\t',header = F,stringsAsFactors = F))

data2 <- list()
for (i in 1:length(data)){
  data2[[i]] <- table(data[[i]]$V1) %>% data.frame()
  data2[[i]]$group <- a[i]
}

df <- do.call('rbind',data2)
df$group %<>% factor(.,levels = c('PQS',paste0('Group',1:6)))

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
my_comparisons = list(c("PQS","Group1"),c("Group1","Group2"),
                      c("Group2","Group3"),c("Group3","Group4"),
                      c("Group4","Group5"),c("Group5","Group6"))
medians <- aggregate(Freq ~  group, df, median)
library(ggpubr)
library(ggplot2)
p <- ggplot(data = df,aes(x = group,y = Freq,fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  #geom_text(data = medians, aes(label = Freq, y = 220),size = 4,col = 'red')+
  scale_x_discrete(limits=c("PQS",paste0('Group',1:6))) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(320,350),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("TF-binding diversity")+
  coord_cartesian(ylim=c(0,400))+
  geom_hline(yintercept = medians[medians$group =='PQS',2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/TF_binding_diversity.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)   

#TF_binding_mean_number
num = c()
for (i in 1:length(data)){
  tmp = (data)[[i]]
  num = c(num,nrow(tmp))
}
df2 <- data.frame(group = a,num = num)
df2$sum <- c(123053,75922,52136,59175,53742,27376,1118855)
df2$mean <- df2$num/df2$sum
df2$group %<>% factor(.,levels = c('PQS',paste0('Group',1:6)))

p <- ggplot(df2 ,aes(x = group, y = mean,fill = group)) +
  geom_bar(stat = 'identity',position = 'stack',width= 0.7,color = 'white') +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,160)) +
  geom_text(aes(label= round(mean),vjust = -0.5), color="black", size=5) +
  cowplot::theme_half_open() +
  ylab('TF-binding mean number')+
  theme(axis.title.x =element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = 'none')

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/TF_binding_mean_number.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)  
# -------------------------------------------------------------------------
sed '1d' m.pqs.G4.primary.txt |awk '$7 != "PQS" {print $0}' > m.G4.primary.txt
/home/qian/source/liftOver ~/G4_analysis/result/savedata2/mouse/m.G4.primary.txt mm10.hg19.rbest.chain.gz mm10Tohg19.G4.bed unmap.bed2
bedtools intersect -a mm10Tohg19.G4.bed -b ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed -wa -wb|awk 'BEGIN{OFS="\t"} {print $4,$11,$14}' > mm10Tohg19.inter.G4.txt

rm(list = ls());gc();rm(list = ls())
df <- read.csv('~/G4_analysis/result/savedata2/TF/mm10Tohg19.inter.G4.txt',sep = '\t',header = F,stringsAsFactors = F) 
df2 <- table(df$V3) %>% data.frame()

group = c(paste0('Group',1:6))
num = c(123053,75922,52136,59175,53742,27376)
df3 <- data.frame(group = group,num = num)

df2$sum <- df3[match(df2$Var1,df3$group),2]
df2$ratio <- df2$Freq/df2$sum

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
ggplot(df2 ,aes(x = group, y = (ratio*100),fill = group)) +
  geom_bar(stat = 'identity',position = 'stack',width= 0.7,color = 'white') +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,2)) +
  geom_text(aes(label= round(ratio*100,2),vjust = -0.5), color="black", size=5) +
  cowplot::theme_half_open() +
  ylab('')+
  theme(axis.title.x =element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = 'none')

# -------------------------------------------------------------------------
awk '$7!=0 {print}' gallus.pqs.G4.txt > gallus.G4.bed
/home/qian/source/liftOver ~/G4_analysis/result/savedata2/gallus/gallus.G4.bed galGal6.hg19.rbest.chain.gz gallusTohg19.G4.bed unmap.bed3
bedtools intersect -a gallusTohg19.G4.bed -b ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed -wa -wb|awk 'BEGIN{OFS="\t"} {print $4,$11,$14}' > gallusTohg19.inter.G4.txt

rm(list = ls());gc();rm(list = ls())
df <- read.csv('~/G4_analysis/result/savedata2/TF/gallusTohg19.inter.G4.txt',sep = '\t',header = F,stringsAsFactors = F) 
df2 <- table(df$V3) %>% data.frame()

group = c(paste0('Group',1:6))
num = c(123053,75922,52136,59175,53742,27376)
df3 <- data.frame(group = group,num = num)

df2$sum <- df3[match(df2$Var1,df3$group),2]
df2$ratio <- df2$Freq/df2$sum

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
ggplot(df2 ,aes(x = group, y = (ratio*100),fill = group)) +
  geom_bar(stat = 'identity',position = 'stack',width= 0.7,color = 'white') +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,8)) +
  geom_text(aes(label= round(ratio*100,2),vjust = -0.5), color="black", size=5) +
  cowplot::theme_half_open() +
  ylab('')+
  theme(axis.title.x =element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = 'none')
# -------------------------------------------------------------------------
awk '$7==0 {print}' gallus.pqs.G4.txt > gallus.pqs.bed
/home/qian/source/liftOver ~/G4_analysis/result/savedata2/gallus/gallus.pqs.bed galGal6.hg19.rbest.chain.gz gallusTohg19.pqs.bed unmap.bed3.1
bedtools intersect -a gallusTohg19.pqs.bed -b ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed -wa -wb |awk 'BEGIN{OFS="\t"} {print $4,$11,"PQS"}' > gallusTohg19.inter.pqs.txt

bedtools intersect -a mm10Tohg19.pqs.bed -b ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed -wa -wb|awk 'BEGIN{OFS="\t"} {print $4,$11,"PQS"}' > mm10Tohg19.inter.pqs.txt
bedtools intersect -a mm10Tohg19.pqs.bed -b gallusTohg19.pqs.bed > mm10Tohg19.inter.gallusTohg19.pqs
bedtools intersect -a ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/pqs.bed -b mm10Tohg19.inter.gallusTohg19.pqs > mm10Tohg19.inter.gallusTohg19.inter.hg19.pqs

bedtools intersect -a mm10Tohg19.G4.bed -b gallusTohg19.G4.bed > mm10Tohg19.inter.gallusTohg19.G4
bedtools intersect -a ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/6groups.bed -b mm10Tohg19.inter.gallusTohg19.G4 > mm10Tohg19.inter.gallusTohg19.inter.hg19.G4

rm(list = ls());gc();rm(list = ls())
mm.pqs <- read.csv('~/G4_analysis/result/savedata2/TF/mm10Tohg19.inter.pqs.txt',sep = '\t',header = F,stringsAsFactors = F) 
mm.g4 <- read.csv('~/G4_analysis/result/savedata2/TF/mm10Tohg19.inter.G4.txt',sep = '\t',header = F,stringsAsFactors = F) 

df <- rbind(mm.g4,mm.pqs)
df2 <- table(df$V3) %>% data.frame()

group = c(paste0('Group',1:6),'PQS')
num = c(123053,75922,52136,59175,53742,27376,1118855)
df3 <- data.frame(group = group,num = num)

df2$sum <- df3[match(df2$Var1,df3$group),2]
df2$ratio <- df2$Freq/df2$sum

pqs <- read.csv('~/G4_analysis/result/savedata2/TF/mm10Tohg19.inter.gallusTohg19.inter.hg19.pqs',sep = '\t',header = F,stringsAsFactors = F) 
g4 <- read.csv('~/G4_analysis/result/savedata2/TF/mm10Tohg19.inter.gallusTohg19.inter.hg19.G4',sep = '\t',header = F,stringsAsFactors = F) 
pqs$V7 <- 'PQS'

df.1 <- rbind(g4,pqs)
df2.1 <- table(df.1$V7) %>% data.frame()
df2.1$sum <- df3[match(df2.1$Var1,df3$group),2]
df2.1$ratio <- df2.1$Freq/df2.1$sum

df2$class <- 'Human & Mouse'
df2.1$class <- 'Human & Mouse & Gallus'
df3 <- rbind(df2,df2.1)

colnames(df3)[1] <- 'group'
df3$group%<>% factor(.,levels = c('PQS',paste0('Group',1:6)))

df3$type <- if_else(df3$group == 'PQS','PQS','G4')

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

ggplot(df3,aes(x = class,y = (ratio*1000),color = group)) +
  geom_line(aes(group = group),size = 1) +
  geom_point(aes(shape = type),size = 4) +
  scale_color_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  coord_cartesian(ylim = c(0,20)) +
  theme(axis.title.x =element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = c(0.8,0.6))
# -------------------------------------------------------------------------
bedtools intersect -a ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4_group_PQS.bed  -b GWAS.hg19.bed -c > G4_group_PQS.GWAS.hg19.txt

rm(list = ls());gc();rm(list = ls())
df <- fread('~/G4_analysis/result/savedata2/TF/G4_group_PQS.GWAS.hg19.txt') %>% data.frame()
#df <- filter(df,V7 %in% c('PQS','G4'))
a = aggregate(df$V8,by = list(df$V7),sum)
num = c(123053,75922,52136,59175,53742,27376,1118855,391404)
df3 <- data.frame(group = c(paste0('Group',1:6),'PQS','G4'),num = num)

a$sum <- df3[match(a$Group.1,df3$group),2]
a$ratio <- a$x/a$sum

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
a$Group.1 %<>% factor(.,levels = c('PQS','G4',c(paste0('Group',1:6)))) 

p <- ggplot(data = a,aes(x = Group.1, y = ratio*100,fill = Group.1)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab(paste0('Average number of SNPs overlapped','\n','with each G4/PQS (x10^3)'))+
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,0.5))+
  geom_text(data = a, aes(label = round(ratio*100,2)),vjust = -1, size = 5)+
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/Average_number_of_SNPs_overlapped.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)  

# -------------------------------------------------------------------------
gunzip cis_eQTLs_all_re.gz
unzip trans_eQTLs_all_re.zip
gunzip trans_eQTLs_all_re.gz

sed '1d' trans_eQTLs_all_re |awk 'BEGIN{OFS="\t"} {print $3,$4,$4+1,$2,$1,$5,$6,$7,"trans"}'|bedtools sort -i > trans.txt
sed '1d' cis_eQTLs_all_re |awk 'BEGIN{OFS="\t"} {print $3,$4,$4+1,$2,$1,$5,$6,$7,"cis"}'|bedtools sort -i > cis.txt
uniq trans.txt > trans.uniq.txt
uniq cis.txt > cis.uniq.txt
cat cis.uniq.txt trans.uniq.txt > cir.trans.txt

scp -r ylxiong@211.69.141.147:/home/ylxiong/G4_database/result/savedata2/eQTL/{cir.txt trans.txt} ./
  
for i in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*bed` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;do j=`basename $i`
bedtools intersect -a $i -b cis.txt -c > cis.${j%.bed};done

for i in `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*bed` /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4.bed;do j=`basename $i`
bedtools intersect -a $i -b trans.txt -c > trans.${j%.bed};done

rm(list = ls());gc();rm(list = ls())
path <- "/home/yulix/G4_analysis/result/savedata2/eQTL"
files <- dir(path)[-c(9,18)]
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

data <- list()
data <- lapply(filepath,function(x)read.csv(x,sep = '\t',header = F,stringsAsFactors = F))

for (i in 1:length(data)){
  data[[i]]$group <- files[i]
}

df <- do.call('rbind',data)
df %<>% separate(group,into = c('regulation','group')) 
df2 <- aggregate(df$V7,by = list(df$group,df$regulation),sum)
df2$Group.1 %<>% gsub('pqs','PQS',.)                   
df2$Group.1 = Hmisc::capitalize(df2$Group.1)

group = c(paste0('Group',1:6),'PQS','G4')
num = c(123053,75922,52136,59175,53742,27376,1118855,391404)
df3 <- data.frame(group = group,num = num)

df2$sum <- df3[match(df2$Group.1,df3$group),2]
df2$ratio <- df2$x/df2$sum

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
df2$Group.1 %<>% factor(.,levels = c('PQS','G4',c(paste0('Group',1:6)))) 

ggplot(data = df2,aes(x = Group.1, y = ratio*100,fill = Group.1)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  #scale_y_continuous(expand = c(0,0))+
  ylab(paste0('Average number of eQTL overlapped','\n','with each G4/PQS (x10^2)'))+
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  #coord_cartesian(ylim = c(0,20))+
  geom_text(data = df2, aes(label = round(ratio*100,2)),vjust = -1, size = 4)+
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)+
  facet_wrap(.~Group.2,scales = 'free')
# -------------------------------------------------------------------------
##GTEx 数据处理合并
scp ylxiong@211.69.141.147:/home/ylxiong/G4_database/result/savedata2/eQTL/GTEx_Analysis_v7_eQTL.tar.gz ./
gzip -d *.gz

rm(list = ls());gc();rm(list = ls())
library(data.table)
path <- "/home/yulix/G4_analysis/result/savedata2/GTEx/GTEx_Analysis_v7_eQTL"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

data <- list()
data <- lapply(filepath,function(x)fread(x))

a = lapply(names(data),function(x)strsplit(x,'.',fixed = T)[[1]][1]) %>% unlist()

list <- lapply(data, function(x)x[,1])
for (i in 1:length(list)){
  list[[i]]$chr <- lapply(list[[i]]$variant_id,function(x)strsplit(x,'_',fixed = T)[[1]][1]) %>% unlist()
  list[[i]]$chr %<>% paste0('chr',.)
  list[[i]]$start <- lapply(list[[i]]$variant_id,function(x)strsplit(x,'_',fixed = T)[[1]][2]) %>% unlist()
  list[[i]]$end <- as.numeric(list[[i]]$start) +1
}

list2 <- lapply(list, function(x)x[,c(2:4,1)])
for (i in 1:length(list)){
  list2[[i]]$name <- a[i]
}
df <- do.call('rbind',list2)
write.table(df,file = '/home/yulix/G4_analysis/result/savedata2/GTEx/GTEx.bed',
            sep = '\t',col.names = F,row.names = F,quote = F)

sed 's/4.1e+07/41000000/g' GTEx.bed |sed 's/1.52e+08/152000000/g' > GTEx.bed2
bedtools sort -i GTEx.bed2 |uniq > GTEx.uniq.bed
bedtools intersect -a ~/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/G4_group_PQS.bed  -b GTEx.uniq.bed -c > G4_group_PQS.GTEx.txt

rm(list = ls());gc();rm(list = ls())
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata2/GTEx/G4_group_PQS.GTEx.txt') %>% data.frame()
df2 <- aggregate(df$V8,by = list(df$V7),sum)

group = c(paste0('Group',1:6),'PQS','G4')
num = c(123053,75922,52136,59175,53742,27376,1118855,391404)
df3 <- data.frame(group = group,num = num)
df2$sum <- df3[match(df2$Group.1,df3$group),2]
df2$ratio <- df2$x/df2$sum

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
df2$Group.1 %<>% factor(.,levels = c('PQS','G4',c(paste0('Group',1:6)))) 

ggplot(data = df2,aes(x = Group.1, y = ratio,fill = Group.1)) +
  geom_bar(stat="identity",width=0.8,position='dodge',color = 'white') +
  scale_y_continuous(expand = c(0,0))+
  ylab(paste0('Average number of GETx overlapped','\n','with each G4/PQS '))+
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none" ) +
  scale_fill_manual(values = color)+
  coord_cartesian(ylim = c(0,1.8))+
  geom_text(data = df2, aes(label = round(ratio,2)),vjust = -1, size = 5)+
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)
# -------------------------------------------------------------------------
#生成 100 bp 的滑动窗口
bedtools makewindows -g genome.chr -w 100 -s 100 > 100bp.bed
bedtools intersect -a /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/100bp.bed -b cis.txt -c > 100bp.inter.cis.txt
bedtools intersect -a /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/100bp.bed -b trans.txt -c > 100bp.inter.trans.txt

#方法二：
#生成 bw 文件
mv 100bp.inter.cis.txt 100bp.inter.cis.bedgraph
/home/qian/source/bedGraphToBigWig 100bp.inter.cis.bedgraph /home/yulix/G4_analysis/data/ref/GRCh37.primary_assembly.genome.fa.fai cis.bw2
mv 100bp.inter.trans.txt 100bp.inter.trans.bedgraph
/home/qian/source/bedGraphToBigWig 100bp.inter.trans.bedgraph /home/yulix/G4_analysis/data/ref/GRCh37.primary_assembly.genome.fa.fai trans.bw2
#compute.pbs
computeMatrix reference-point -S cis.bw2 -R `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` -b 100 -a 100 -o cis.mat.gz2
computeMatrix reference-point -S trans.bw2 -R `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` -b 100 -a 100 -o trans.mat.gz2

computeMatrix reference-point -S cis.bw2 -R `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` -b 1000 -a 1000 -o cis.mat.gz3
computeMatrix reference-point -S trans.bw2 -R `ls /home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group/*` -b 1000 -a 1000 -o trans.mat.gz3
#画图
plotProfile -m cis.mat.gz3  -out cis.3.pdf --colors "#FDDBC7" "#DFAFA5" "#C18383" "#A35762" "#842B40" "#67001F" "#4d4d4d"
plotProfile -m trans.mat.gz3  -out trans.3.pdf --colors "#FDDBC7" "#DFAFA5" "#C18383" "#A35762" "#842B40" "#67001F" "#4d4d4d"
#脚本
cd /home/ylxiong/G4_database/result/savedata2/eQTL
computeMatrix reference-point -S cis.bw2 -R `ls ../group/*` -b 10000 -a 10000 -o cis.mat.gz
plotProfile -m cis.mat.gz  -out cis.10K.pdf --colors "#FDDBC7" "#DFAFA5" "#C18383" "#A35762" "#842B40" "#67001F" "#4d4d4d"

cd /home/ylxiong/G4_database/result/savedata2/eQTL
computeMatrix reference-point -S trans.bw2 -R `ls ../group/*` -b 10000 -a 10000 -o trans.mat.gz
plotProfile -m trans.mat.gz  -out trans.10K.pdf --colors "#FDDBC7" "#DFAFA5" "#C18383" "#A35762" "#842B40" "#67001F" "#4d4d4d"
qsub -q batch -V -l nodes=1:ppn=4 trans.pbs
qsub -q batch -V -l nodes=1:ppn=4 cis.pbs

# G4_correlation-------------------------------------------------------------------------
cat list| while read id;do arr=($id);a=${arr[0]};b=${arr[1]};bedtools reldist -a ../pqs.intersect.narrowpeak/group/$a -b ../pqs.intersect.narrowpeak/group/$b > ${a%.bed}_${b%.bed}.reldist;done
for i in *reldist;do echo $i;done >1
for i in *reldist;do awk '{print $1*$4}' $i |awk '{sum += $1};END {print sum}';done >> res
paste res 1 -d "\t" |sed s/.reldist//g > result

rm(list = ls());gc();rm(list = ls())
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata2/G4_correlation/result')
df <- separate(df,col = V2,sep = '_',into = c('x','y'))
df$V1 <- 1-(df$V1)
df <- spread(df,key = y,value = V1)
df <- column_to_rownames(df,var = 'x')
df[is.na(df)] <- 1
library(Hmisc)
colnames(df) <- capitalize(colnames(df))
rownames(df) <- capitalize(rownames(df))

a <- data.frame()
for (i in 1:6){
  for (j in 1:6){
    a[i,j] <- (df[i,j] + df[j,i])/2
  }
}
colnames(a) <- c(paste0('Group',1:6))
rownames(a) <- c(paste0('Group',1:6))
library(corrplot)
# 混合方法之上三角为圆形，下三角为数字
col <- colorRampPalette(colors =c('#2166ac','white','#b2182b'))
a = as.matrix(a)
corrplot(a,type="upper",tl.pos="tp",tl.col = 'black',col = col(20))
corrplot(a,add=TRUE, type="lower",
         method="number",diag=FALSE,tl.pos="n", cl.pos="n",tl.col = 'black',col = col(20))
 
# group富集的TF之间的相关性-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
library(tidyverse)
path <- '/home/yulix/G4_analysis/result/savedata2/TF/ChIP-Atlas'
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})                                                  
data <- lapply(filepath,function(x)fread(x))       

data2 <- lapply(data,function(x)x[,c(3,9,11)])  

data3 <- lapply(data2,function(x)x[x$V11 != 'Inf',])  

data4 <- lapply(data3,function(x)x[x$V9 < -5 & x$V11 > 3 ,])  
df <- list()
for (i in 1:6){
  df[[i]] <- tapply(data4[[i]]$V11,data4[[i]]$V3,mean) 
}

names(df) <- names(data4)

df2 <- list()
for (i in 1:6){
  df2[[i]] <- df[[i]] %>% data.frame() %>% rownames_to_column(.)
}

names(df2) <- names(df)

group1 <- df2[[1]]
group2 <- df2[[2]]
group3 <- df2[[3]]
group4 <- df2[[4]]
group5 <- df2[[5]]
group6 <- df2[[6]]

df3 <- do.call("rbind",df2) %>% data.frame()
df3 <- data.frame(df3$rowname) %>% unique()
colnames(df3)[1] <- 'antigen_name'

df3$group1 <- group1[match(df3$antigen_name,group1$rowname),2]
df3$group2 <- group2[match(df3$antigen_name,group2$rowname),2]
df3$group3 <- group3[match(df3$antigen_name,group3$rowname),2]
df3$group4 <- group4[match(df3$antigen_name,group4$rowname),2]
df3$group5 <- group5[match(df3$antigen_name,group5$rowname),2]
df3$group6 <- group6[match(df3$antigen_name,group6$rowname),2]

df3[is.na(df3)] <- 0
rownames(df3) <- df3$antigen_name
df4 <- df3[,2:7]
library(Hmisc)
colnames(df4) <- capitalize(colnames(df4))
rownames(df4) <- capitalize(rownames(df4))
# 相关系数图
library(corrplot)
col <- colorRampPalette(c('#2166ac','white','#b2182b'))
corrplot(corr =cor(df4[1:6]),type="upper",tl.pos="tp",tl.col = 'black',col = col(20))
corrplot(corr = cor(df4[1:6]),add=TRUE, type="lower",
         method="number",diag=FALSE,tl.pos="n", cl.pos="n",tl.col = 'black',col = col(20))

# 画富集前10 的 TF 表格 -------------------------------------------------------------------------
df2.1 = lapply(df2, function(x)x[order(-x$.),][1:10,1])
df2.2 <- do.call('cbind',df2.1) %>% data.frame()
colnames(df2.2) <- c(paste0('Group',1:6))

library(ggpubr)
tbody.style = tbody_style(color = "black",
                          fill = c("#e8f3de", "#d3e8bb"), hjust=1, x=0.9)
p <- ggtexttable(df2.2, rows = NULL,
                 theme = ttheme(
                   colnames.style = colnames_style(color = "white", fill = "#8cc257"),
                   tbody.style = tbody.style))

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/TF_table.pdf",
                  width = 5,height = 3.5,devices = NULL,res = 300,show = F)    
# cell_type-------------------------------------------------------------------------
ls /home/ylxiong/G4_database/newdata/Homo_sapiens/CHIP-seq/* |grep -v input |grep -v ^[l,s,d] |grep -v / |sed '/^$/d' |cut -d "_" -f1 |sort |uniq -c > newdata.num
ls /home/ylxiong/G4_database/data/CHIP-seq/*/*/ |grep -v ^SRR |grep -v input > data.txt

rm(list = ls());gc();rm(list = ls())
df <- fread('/home/yulix/G4_analysis/result/savedata2/cell_type/data.txt') %>% data.frame()
df2 <- fread('/home/yulix/G4_analysis/result/savedata2/cell_type/newdata.num') %>% data.frame()
colnames(df2) <- c('V2','V1')

df3 <- rbind(df,df2)
df4 <- aggregate(df3$V2,by = list(df3$V1),sum) %>% data.frame()
df4[14,2] <- 12
#write.table(df4,file = '/home/yulix/G4_analysis/result/savedata2/cell_type/sample.txt',
            sep = '\t',col.names = F,row.names = F,quote = F)

a <- df4[order(df4$x),1]
df4$Group.1 %<>% factor(.,levels = a)
p <- ggplot(df4 ,aes(x = Group.1, y = x)) +
  geom_bar(stat = 'identity',position=position_dodge(1),width=1,color = 'white',fill ='#d0b8b9')+
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=12),
        axis.title.x = element_text(size = 16),
        axis.title.y =element_blank(),legend.position="none")+
  ylab("Number of samples")+
  coord_flip(ylim = c(0,30))  
Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/cell_type.pdf",
                  width = 6.8,height = 6.8,devices = NULL,res = 300,show = F) 
91bfdb
# G4_peak_density-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata2/cell_type/G4peak.length') %>% data.frame()

library(ggplot2)
n = c(nrow(df),median(df$V1))

p <- ggplot(df,aes(x=log10(V1+1)))+geom_density(fill="#D7AF9E")+
  cowplot::theme_half_open()+
  theme(axis.title = element_text(size = 16),
        axis.text= element_text(size = 14))+
  ylab("Density (%)")+
  xlab(bquote(log[10](length+1)))+ 
  annotate(geom="text", x=3, y=1, size=4,
           label=paste0("Total G4 Peaks","\n",n[1],"\n","Median length","\n",n[2],' bp'))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  coord_cartesian(xlim = c(1.5,3.5))

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/G4_peak_density.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# G4 peak 数量分布图-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
num <- fread("/home/yulix/G4_analysis/result/savedata2/cell_type/narrowPeak.num")
colnames(num) <- c("num","sample")   
library(ggplot2)
p<- ggplot(data = num,aes(x = sample, y = num,fill = sample)) +
  geom_bar(stat="identity",width=1,position='dodge',fill = '#D7AF9E',color ='white') +
  scale_y_continuous(expand = c(0,0))+
  cowplot::theme_half_open() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  ylab('Number') +
  xlab('Sample')

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/G4_peak_number.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

#  ATAC peak 数量分布图-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
num <- fread("/home/yulix/G4_analysis/result/savedata2/cell_type/atac.num")
colnames(num) <- c("num","sample")   
library(ggplot2)
p2 <- ggplot(data = num,aes(x = sample, y = num,fill = sample)) +
  geom_bar(stat="identity",width=1,position='dodge',fill = '#D7AF9E',color ='white') +
  scale_y_continuous(expand = c(0,0))+
  cowplot::theme_half_open() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  ylab('Number') +
  xlab('Sample') + 
  coord_cartesian(ylim = c(0,150000))
Ipaper::write_fig(p2,file = "/home/yulix/G4_analysis/result/figure2/ATAC_peak_number.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# -------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
gene.exp <- fread("/home/yulix/G4_analysis/result/savedata2/RNA_seq/gene_exp.bed") %>% data.frame()
gene.exp$mean <- rowMeans(gene.exp[,5:45])
df <- gene.exp[,c(1:4,46)]

G4.num <- fread('/home/yulix/G4_analysis/result/savedata/RNA_seq/SRR11113121_G4') %>% data.frame()
G4.length <- fread('/home/yulix/G4_analysis/result/savedata/RNA_seq/SRR11113121_G4_length') %>% data.frame()

df$num <- G4.num[match(df$gene_id,G4.num$V4),6]

a <- aggregate(G4.length$V13,by = list(G4.length$V4,G4.length$V12),sum) %>% data.frame()
colnames(a) <- c('gene_id','group','length')
length <- tapply(a$length, a$gene_id, sum) %>% data.frame()
length <- rownames_to_column(length)
colnames(length) <- c('gene_id','length')
df$length <- length[match(df$gene_id,length$gene_id),2]
df[is.na(df)] <- 0
df$gene_length <- df$end - df$start
df$density <- df$length/df$gene_length 

cor.test(df$density,df$mean)
cor.test(df$length,df$mean)
cor.test(df$num,df$mean)  

p1 <- ggplot(data = df,aes(x = length,y = mean)) +
  geom_point(shape = 21,size = 1.5,color = '#D7AF9E',fill = 'white',alpha = 1) +
  cowplot::theme_half_open()+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  xlab("Length")+
  ylab("Expression level (TPM)")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))+
  coord_cartesian(ylim = c(0,3000),xlim = c(0,4000)) +
  geom_smooth(method = 'gam',color = 'black') +
  annotate("text",x = 1000,y = 2000,label = "r = 0.05;pvalue = 2.2e-16",size = 5)

Ipaper::write_fig(p1,file = "/home/yulix/G4_analysis/result/figure2/length_expression_poinplot.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)     
p2 <- ggplot(data = df,aes(x = num,y = mean)) +
  geom_point(shape = 21,size = 1.5,color = '#D7AF9E',fill = 'white',alpha = 1) +
  cowplot::theme_half_open()+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  xlab("Number")+
  ylab("Expression level (TPM)")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))+
  coord_cartesian(ylim = c(0,5000),xlim = c(0,500)) +
  geom_smooth(method = 'gam',color = 'black') +
  annotate("text",x = 200,y = 3000,label = "r = 0.05;pvalue = 2.2e-16",size = 5)
Ipaper::write_fig(p2,file = "/home/yulix/G4_analysis/result/figure2/num_expression_poinplot.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)  
p3 <- ggplot(data = df,aes(x = density,y = mean)) +
  geom_point(shape = 21,size = 1.5,color = '#D7AF9E',fill = 'white',alpha = 1) +
  cowplot::theme_half_open()+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  xlab("Density")+
  ylab("Expression level (TPM)")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))+
  coord_cartesian(ylim = c(0,3000),xlim = c(0,0.4)) +
  geom_smooth(method = 'gam',color = 'black') +
  annotate("text",x = 0.2,y = 2000,label = "r = 0.06;pvalue = 2.2e-16",size = 5)
Ipaper::write_fig(p3,file = "/home/yulix/G4_analysis/result/figure2/density_expression_poinplot.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)  

# gallusTohg19-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
pqs <- read.csv('~/G4_analysis/result/savedata2/TF/gallusTohg19.inter.pqs.txt',sep = '\t',header = F,stringsAsFactors = F) 
g4 <- read.csv('~/G4_analysis/result/savedata2/TF/gallusTohg19.inter.G4.txt',sep = '\t',header = F,stringsAsFactors = F) 

df <- rbind(g4,pqs)
df2 <- table(df$V3) %>% data.frame()

group = c(paste0('Group',1:6),'PQS')
num = c(123053,75922,52136,59175,53742,27376,1118855)
df3 <- data.frame(group = group,num = num)

df2$sum <- df3[match(df2$Var1,df3$group),2]
df2$ratio <- df2$Freq/df2$sum

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

df2$Var1 %<>% gsub('PQS','non-eG4',.)
df2$Var1 %<>% factor(.,levels = c("non-eG4",paste0('Group',1:6)))

p <- ggplot(df2 ,aes(x = Var1, y = (ratio*100),fill = Var1)) +
  geom_bar(stat = 'identity',position = 'stack',width= 0.7,color = 'white') +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,8)) +
  geom_text(aes(label= round(ratio*100,2),vjust = -0.5), color="black", size=5) +
  cowplot::theme_half_open() +
  ylab('(%) Conserved at structure level')+
  theme(axis.title.x =element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = 'none')

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/gallusTohg19.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F) 

# length.non-eG4_group-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/group'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

a = c(paste0('Group',1:6),'Non-eG4')
for (i in 1:length(data)){
  data[[i]]$group <- a[i]
}

df <- do.call('rbind',data)
df$length <- df$V3-df$V2
medians <- aggregate(length ~  group, df, median)
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
my_comparisons = list(c("Non-eG4","Group1"),c("Group1","Group2"),
                      c("Group2","Group3"),c("Group3","Group4"),
                      c("Group4","Group5"),c("Group5","Group6"))

library(ggpubr)
df$group %<>% factor(.,levels = c("Non-eG4",paste0('Group',1:6)))

p <- ggplot(data = df,aes(x = group,y = length,fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  #geom_text(data = medians, aes(label = length, y = 40),size = 4,col = 'red')+
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(60,70),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("Length (bp)")+
  coord_cartesian(ylim=c(0,80))+
  geom_hline(yintercept = medians[medians$group =='Non-eG4',2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/length.non-eG4_group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)
# GC_content.non-eG4_group-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
path <- '/home/yulix/G4_analysis/result/savedata2/GC_content'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

a = c(paste0('Group',1:6),'Non-eG4')
for (i in 1:length(data)){
  data[[i]]$group <- a[i]
}

df <- do.call('rbind',data)
colnames(df)[6] <- 'gc'
medians <- aggregate(gc ~  group, df, median)

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
my_comparisons = list(c("Non-eG4","Group1"),c("Group1","Group2"),
                      c("Group2","Group3"),c("Group3","Group4"),
                      c("Group4","Group5"),c("Group5","Group6"))

library(ggpubr)
df$group %<>% factor(.,levels = c("Non-eG4",paste0('Group',1:6)))

p <- ggplot(data = df,aes(x = group,y = gc,fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  #geom_text(data = medians, aes(label = round(gc,3), y = 0.5),size = 4,col = 'red')+
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(1,1.1),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("GC content")+
  coord_cartesian(ylim=c(0.4,1.2))+
  geom_hline(yintercept = medians[medians$group =='Non-eG4',2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/GC_content.non-eG4_group.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)
# fraction.G4_group-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
df <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/fraction.G4_group.txt') %>% data.frame()
df$Class %<>% gsub('G4-depleted genes','eG4-depleted genes',.) %>% gsub('G4-containing genes','eG4-containing genes',.)

color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)
df$Class %<>% factor(.,levels = c('eG4-depleted genes','eG4-containing genes',c(paste0('Group',1:6))))  

G4.num <- fread('/home/yulix/G4_analysis/result/savedata2/RNA_seq/G4.num.txt') %>% data.frame()
#nrow(G4.num[G4.num$class == 'G4-depleted genes' & G4.num$expression > 1,])/table(G4.num$class)[[2]]  #0.1974195
#nrow(G4.num[G4.num$class == 'G4-containing genes' & G4.num$expression > 1,])/table(G4.num$class)[[1]] #0.532067
a = data.frame(x = c(nrow(G4.num[G4.num$class == 'G4-depleted genes' & G4.num$expression > 1,]),table(G4.num$class)[[2]]),
               y = c(nrow(G4.num[G4.num$class == 'G4-containing genes' & G4.num$expression > 1,]),table(G4.num$class)[[1]]))
fisher.test(a)

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
  geom_vline(xintercept = 2.5,linetype = 2,size=0.8)+
  annotate("text",x = 1.5,y = 0.8,label = "p < 2.2e-16",size = 4)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/fraction.G4_group.pdf",
                  width = 7.5,height = 6.5,devices = NULL,res = 300,show = F)      

# fold_enrichment.desease-------------------------------------------------------------------------
df <- fread('/home/yulix/G4_analysis/result/savedata/G4_TE/disease/fold_enrichment.desease.txt') %>% data.frame()

df$class %<>% factor(.,levels = c('Lethal','Lethal_disease','Disease','OtherGenes','Oncogene','TSG'))
df$group %<>% gsub('G4','eG4',.) %>% gsub('PQS','Non-eG4',.)
df$group %<>% factor(.,levels = c('Non-eG4','eG4',paste0('Group',1:6)))
color = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c(c('#4d4d4d','#b2182b'),color)

#计算显著性
a_obs <- gather(sum_group_obs,key = type,value = sum_group_obs,2)
a_simu <- gather(sum_group_simu,key = type,value = sum_group_simu,2)
colnames(a_obs)[4] <- 'value'
colnames(a_simu)[4] <- 'value'
a = rbind(a_obs,a_simu)
a <- filter(a,group %in% c('G4','PQS'))
a1 <-spread(data=a, key=group,value=value)
Disease <- a1[a1$class == 'Disease',3:4]
fisher.test(Disease) #p-value < 2.2e-16

Lethal <- a1[a1$class == 'Lethal',3:4]
fisher.test(Lethal) #p-value < 2.2e-16

Lethal_disease <- a1[a1$class == 'Lethal_disease',3:4]
fisher.test(Lethal_disease) #p-value < 2.2e-16

Oncogene <- a1[a1$class == 'Oncogene',3:4]
fisher.test(Oncogene) #p-value < 2.2e-16

OtherGenes  <- a1[a1$class == 'OtherGenes',3:4]
fisher.test(OtherGenes ) #p-value < 2.2e-16

Tsg <- a1[a1$class == 'Tsg',3:4]
fisher.test(Tsg) #p-value < 2.2e-16

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
  geom_hline(yintercept = 1,linetype = 2) +
  annotate("text",x = c(0.7,1.7,2.7,3.7,4.8,5.7),y = c(1.6,1.6,1.6,1.7,2.1,1.5),
           label = "p < 2.2e-16",size = 3)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/fold_enrichment.desease.pdf",
                  width = 7,height = 4,devices = NULL,res = 300,show = F)   

# Reproduction_rate-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(data.table)
df <- fread('/home/yulix/G4_analysis/result/savedata2/pqs.intersect.narrowpeak/hg.pqs.result') %>% data.frame()

U2OS <- df[,c("U2OS_hypoxia","U2OS_normoxia")]
U2OS.1 <-  mutate(U2OS,sum = rowSums(U2OS[,1:2]))
table(U2OS.1$sum)[[3]]/(table(U2OS.1$sum)[[3]] + table(U2OS.1$sum)[[2]]) #0.4767506

hacat <- df[,c("HACAT_NA","HACAT_Entinostat")]
hacat.1 <-  mutate(hacat,sum = rowSums(hacat[,1:2]))
table(hacat.1$sum)[[3]]/(table(hacat.1$sum)[[3]] + table(hacat.1$sum)[[2]]) #0.4827727

k562 <- df[,c("K562_DRB","K562_hypoxia","K562_NA","K562_noDRB","K562_normoxia","K562_noTPL","K562_TPL")]

x = c(rep(1,6),rep(2,5),rep(3,4),rep(4,3),rep(5,2),rep(6,1))
y = c(2:7,3:7,4:7,5:7,6:7,7)
a_k562 = data.frame(x = x,y = y)

b_k562 = vector()
for (i in 1:21){
  k562.1 <- k562[,c(a_k562[i,1],a_k562[i,2])]
  k562.1 <-  mutate(k562.1,sum = rowSums(k562.1[,1:2]))
  b_k562[i] = table(k562.1$sum)[[3]]/(table(k562.1$sum)[[3]] + table(k562.1$sum)[[2]])}

HEK293T <- df[,c("HEK293T_CUT.Tag","HEK293T_Flavopiridol.CUT.Tag_DMSO",
                 "HEK293T_Flavopiridol.CUT.Tag_FP","HEK293T_PDS.CUT.Tag_DMSO",         
                 "HEK293T_PDS.CUT.Tag_PDS","HEK293T_TMPyP4.CUT.Tag_DMSO",      
                 "HEK293T_TMPyP4.CUT.Tag_TMPYP4")]

b_HEK293T = vector()
for (i in 1:21){
  HEK293T.1 <- HEK293T[,c(a_k562[i,1],a_k562[i,2])]
  HEK293T.1 <-  mutate(HEK293T.1,sum = rowSums(HEK293T.1[,1:2]))
  b_HEK293T[i] = table(HEK293T.1$sum)[[3]]/(table(HEK293T.1$sum)[[3]] + table(HEK293T.1$sum)[[2]])}

within = c(0.4767506,0.4827727,b_k562,b_HEK293T)

a = data.frame(U2OS,k562)
b = data.frame(hacat,HEK293T)

c = list()
for (i in 1:100){
  set.seed(i)
  c[[i]] = sample(9,2,replace = T)}

c.df <- do.call('rbind',c) %>% data.frame()
c.df <- unique(c.df)

between = vector()
for (i in 1:nrow(c.df)){
  ab <- data.frame(a[,c.df[i,1]],b[,c.df[i,2]])
  ab <-  mutate(ab,sum = rowSums(ab[,1:2]))
  between[i] = table(ab$sum)[[3]]/(table(ab$sum)[[3]] + table(ab$sum)[[2]])}

within <- data.frame(class = 'within',ratio = within)
between <- data.frame(class = 'between',ratio = between)

df2 <- rbind(within,between)
df2$class <- Hmisc::capitalize(df2$class)
df2$class %<>% factor(.,levels = c("Within","Between"))

medians <- aggregate(ratio ~  class, df2, median)

p <- ggplot(data = df2,aes(x = class,y = ratio,fill = class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white",alpha = 0.8) +
  #geom_text(data = medians, aes(label = round(score,3), y = 85),size = 8,col = 'red')+
  ggpubr::stat_compare_means(comparisons = list(c("Within","Between")),label.y = 0.8,
                             aes(label = paste0("p = ", ..p.format..)),
                             tip.length = 0,size=5)+
  scale_fill_manual(values = c('#e9a3c9','#bababa') ) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  ylab("Reproduction rate of G4 in cells")+
  coord_cartesian(ylim=c(0,1))+
  geom_hline(yintercept = medians[1,2],linetype=2)

Ipaper::write_fig(p,file = "/home/yulix/G4_analysis/result/figure2/Reproduction_rate.pdf",
                  width = 6,height = 4,devices = NULL,res = 300,show = F)  
































