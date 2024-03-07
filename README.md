# TFenrichment_Yu

Source codes used for the ChIP-seq-based transcription factor enrichment analysis.

Instruction
```
library(stringr)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(dplyr)
library(tibble)
```

The file corresponding to the Transcription Factor Target Matrix and TF lists are stored in AnnoList2.RData and peakAnnoNames.RData.


The TF list was compiled based on the referenced paper by 
> Lambert, S. A., Jolma, A., Campitelli, L. F., Das, P. K., Yin, Y., Albu, M., Chen, X., Taipale, J., Hughes, T. R., & Weirauch, M. T. (2018). The Human Transcription Factors. Cell, 172(4), 650–665. https://doi.org/10.1016/j.cell.2018.01.029

The Lysosomal_biogenesis_Bordi gene list was compiled based on the referenced paper by
> Bordi, M., De Cegli, R., Testa, B. et al. A gene toolbox for monitoring autophagy transcription. Cell Death Dis 12, 1044 (2021). https://doi.org/10.1038/s41419-021-04121-9

```
System info:
R version 4.2.3 (2023-03-15 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.1.2                              ChIPseeker_1.34.1                        TxDb.Hsapiens.UCSC.hg38.knownGene_3.16.0 GenomicFeatures_1.50.4                  
 [5] AnnotationDbi_1.60.2                     Biobase_2.58.0                           GenomicRanges_1.50.2                     GenomeInfoDb_1.34.9                     
 [9] IRanges_2.32.0                           S4Vectors_0.36.2                         BiocGenerics_0.44.0                      tidyr_1.3.0                             
[13] stringr_1.5.0                           
```

no installation is needed

demo files in sample_bed folder

Expected run time for demo :~1min
