---
title: "rna-seq-volcano_plots"
output: 
  html_document:
    keep_md: true
---



Loading libraries


Read in expression data

```r
path=readtext("C:/Users/jwang/Desktop/DiffExpression-path.txt") # get filename from path


df <- read_excel(path$text, sheet = 2)
df <- arrange(df, desc(abs(stat)))#noduplicates
df <- df[!duplicated(df$GeneID),]
#add boundaries for logfc and pval
df$fdr<- p.adjust(df$pvalue, method="fdr")

df <- dplyr::rename(df, logFC= colnames(df)[4])

#set boundaries on logfc and pvalue
pthresh= 10^-8
fcthresh= 6

df$logFC[df$logFC>fcthresh]= fcthresh
df$logFC[df$logFC< -fcthresh]= -fcthresh
df$pvalue[df$pvalue< pthresh]= pthresh
head(df[,1:9])
```

```
## # A tibble: 6 x 9
##   GeneID     GeneName IsCoding  logFC lfcSE  stat   pvalue     padj healthy_s110
##   <chr>      <chr>    <chr>     <dbl> <dbl> <dbl>    <dbl>    <dbl>        <dbl>
## 1 ENSG00000~ ERAP2    TRUE      5.24  0.476 11.0   1.00e-8 5.99e-24         20.6
## 2 ENSG00000~ LAMA1    TRUE      2.75  0.294  9.36  1.00e-8 8.01e-17         57.7
## 3 ENSG00000~ IGFBP4   TRUE      1.06  0.170  6.24  1.00e-8 2.51e- 6       7064. 
## 4 ENSG00000~ ANKRD2   TRUE     -3.03  0.487 -6.22  1.00e-8 2.51e- 6        657. 
## 5 ENSG00000~ HOXC10   TRUE      0.921 0.166  5.56  2.73e-8 1.07e- 4        634. 
## 6 ENSG00000~ KLHL4    TRUE     -0.732 0.135 -5.41  6.28e-8 1.93e- 4        786.
```


```r
p <- EnhancedVolcano(df, lab=rep("",nrow(df)), x= "logFC", y="fdr",pCutoff=.05, FCcutoff=1, pointSize = 1.5, labSize = 3.0, xlab = bquote(~Log[2]~ "Fold Change (IO/Healthy)"), ylab = bquote(~-Log[10]~italic("FDR")), ylim=c(0,5), legend = c("d","d","d","FDR<.05 |FC|>1"), title="",subtitle="",col=c('black','red','royalblue','green')          ) #,colCustom = keyvals  #for mouse coloring

p+ ggplot2::coord_cartesian(xlim=c(-6, 6)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-6,6, 2))

```

![volcano-plot](https://user-images.githubusercontent.com/46359281/74269596-b7969480-4cd7-11ea-9185-31b58d334bdf.png)

![](volcano_plots-rna-seq_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


