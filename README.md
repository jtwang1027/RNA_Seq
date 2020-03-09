# RNA Seq
Here we start with an expression matrix generated from DESeq2 and perform various analyses (gene set enrichment analysis and correlation tests) as well as generate common visualizations (heatmaps, venn diagrams, volcano plots, violin plots).


# *Generating heatmaps*
Starting from the expression matrix (ie output from DESeq2 or EdgeR):

1. Select genes to show: (for example using FC and pvalue cutoff or lfcShrink cutoff)
2. Z-scaling each gene
3. Calling heatmap.2 

![mouse-heatmap-3groups_reversal-alt](https://user-images.githubusercontent.com/46359281/76250075-f4907100-621a-11ea-956a-43688a2bfb2d.jpg)

    

# *Generating Volcano plots using EnhancedVolcano*
Source code [(link)](https://github.com/jtwang1027/RNA_Seq/blob/master/volcano_plots-rna-seq.md) 
![volcano-plot](https://user-images.githubusercontent.com/46359281/74269596-b7969480-4cd7-11ea-9185-31b58d334bdf.png)

  
  

# *Generating Venn Diagrams*

![mouse_treatment_venn_fdr05_logfc1-testest](https://user-images.githubusercontent.com/46359281/76251573-9913b280-621d-11ea-9ece-aacc2d2bc4cf.png)

