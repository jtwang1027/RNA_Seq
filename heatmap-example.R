#https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
#https://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
#https://blog.bioturing.com/2018/09/24/heatmap-color-scale/
#core used z-score plotting, not logFC
require(dplyr)
# require(fgsea)
# require(biomaRt)
require(gplots)
require(apeglm)
require(readxl)
require(RColorBrewer)
require(stringr)

###################REVERSED GENES FROM therapeutic TREATMENT
#30 genes (fdr <.05), generated from: 'ms_venn_diagram.rmd'
overlap=read.csv('C:/Users/jwang/Documents/DUKE/LAB/180710-wt-aav-ko-rna-seq-koeberl-RNA-SEQUENCINGmouse-aav_reversed-fdr05.txt',sep='\t', header= FALSE, stringsAsFactors = FALSE)
overlap=overlap[[1]]
###################

genes_shown <-15

#first dataset
genes <- read_excel("C:/path/to/file.xlsx",   sheet='KO.vs.WT') %>% as.data.frame()
genes=genes[ , !(names(genes) %in% "logFC (KO / WT)")] #critical for how columns are selected for heatmap later
genes=genes[genes$GeneName %in% overlap,]

#2nd dataset
genes2 <- read_excel("path/to/file2.xlsx",sheet='AAV.vs.KO') %>% as.data.frame() 
genes2=genes2[genes2$GeneName %in% overlap,c('stat',"GeneName","AAV.G5.315","AAV.G5.316", "AAV.G5.317","AAV.G5.318")]

mer=merge(genes, genes2, by= 'GeneName')
mer=mer[mer$IsCoding=='TRUE',]

mer=mutate(mer, stat_comb=abs(stat.x)+abs(stat.y)) #combined stat score
mer <- mer %>% arrange(desc(stat_comb))
mer= filter(mer, IsCoding=='TRUE')


#split to top and bottom
idx_pos=which(mer$stat.x>0)[1:(genes_shown/2)]
idx_neg=which(mer$stat.x<0)[1:(genes_shown/2)]

heat_mat2 <- mer[unique(c(idx_pos ,idx_neg, valid_idx)), str_detect(colnames(mer), 'WT|KO|AAV')]
row.names(heat_mat2) <-mer[unique(c(idx_pos ,idx_neg, valid_idx)),'GeneName'] # needed for heatmap later

#reorient for it will be WT on top when rotated
heat_mat2 <- heat_mat2[,c(9:12,5:8,1:4)]
colnames(heat_mat2)= c( paste0("AAV",1:4), paste0("KO",1:4),paste0("WT",1:4) )


heat_mat2 <- t(scale(t(heat_mat2), center = T, scale = T)) # transpose and normalize, Z-score from all genes
#you need to reverse the normal red-->blue gradient, so that red is up
myPalette <- colorRampPalette(rev(brewer.pal(11,"RdBu")))

{
  dev.off()
  heatmap.2(heat_mat2, 
            col= myPalette(100),
            Rowv=T,
            Colv=F,
            trace="none", 
            lwid= c(1.5,1),
            lhei= c(5,8),
            dendrogram = "row",
            density.info="none",
            cexRow = .6,
            cexCol = .8,
            scale = "none",
            key=T,
            lmat=rbind( c(4,3),c(2,1)),
            sepwidth=c(0.02,0.02),
            sepcolor="black",
            rowsep=0:nrow(heat_mat2),
            colsep=0:ncol(heat_mat2),
            margins=c(3,6),
            offsetRow = -.5,
            offsetCol = 1,
            srtRow = -45
  )
  
  
  
}




