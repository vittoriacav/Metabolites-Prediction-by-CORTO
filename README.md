# Predicting metabolites using networks and gene expression profiles from the Cancer Cell Line Encyclopedia

As one of the most recent -omics sciences, metabolomics revealed to have important implications in the oncology field. However, metabolomic studies still have many limits related to the cost of its technologies and the heterogeneity of metabolites that are measured in different studies. One possibility to overcome these problems is to predict metabolites abundances starting from transcript levels measured in the same samples. In fact, compared to metabolomics, transcriptomics is more established, and the related technologies are far more accurate and cheaper. In this thesis project an artificial intelligence algorithm, called *corto*, has been used to test the possibility of metabolites prediction, starting from co-measured transcript abundances. This kind of analyses needs large-scale datasets and, in this case, the recent Cancer Cell Line Encyclopedia (CCLE) dataset, which co-measure hundreds of metabolites and thousands of transcripts, from more than 20 cancer types, was employed. In addition, four more datasets have been collected and have been employed as test sets to confirm predictions based on the CCLE relationships. Among the four additional datasets, two were based on cell lines while the other on patient samples. The results showed a significant capability to predict metabolites from transcript data, with significant correlations between predicted and measured metabolites for more than 20 molecular species. We could predict metabolic levels even in patient datasets for some molecules, such as the uracil, but the best performance was obtained when using cell lines data. Compared to other predicting tools, *corto* showed better or equal performance in prediction but less memory usage and computational time. We believe that *corto*, while not able to fully replace future metabolomics measurements, can be used as a formidable tool to predict metabolite abundances in contexts where only transcript data is present, due to technical reasons or budget limitations. 
![image](https://user-images.githubusercontent.com/76695084/156879046-2c4d89a6-2c2b-48fa-a800-15832b5fb864.png)
