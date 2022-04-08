# Task 2: Calculate and interpret similarity results

There are 4 narrowPeak files stored under /data. These are ChIP-seq datasets, for two cell types (HepG2 and HeLaS3) and two transcription factors (Jun and CTCF). I filtered the data to keep things smaller and faster for this assignment.

Use your functions to calculate the pairwise Jaccard similarity for this set of 4 experiments, then answer the following questions:

pairwise Jaccard similarity:

|               |helas3_ctcf|helas3_jun|hepg2_ctcf|hepg2_jun|
|---------------|-----------|----------|----------|---------|
|**helas3_ctcf**|   1.000   |  0.020   | 0.610    | 0.023   |
|**helas3_jun** |   0.020   |  1.000   | 0.013    | 0.170   |
|**hepg2_ctcf** |   0.610   |  0.013   | 1.000    | 0.026   |
|**hepg2_jun**  |   0.023   |  0.170   | 0.026    | 1.000   |

1. Which two interval sets are the most similar? (10 pts) 
**CTCF ChIP-seq of HepG2 and HeLaS3 are the most similar interval sets (other than themselves) with Jaccard similarity of 0.610.**

2. Which two interval sets are the most different? (10 pts) 
**CTCF ChIP-seq of HepG2 and Jun ChIP-seq of HeLaS3 are the most different interval sets with Jaccard similarity of 0.013.** 

3. Based on these results, which factor, CTCF or Jun, would you predict varies more across cell types? (10 pts)
**Based on the results, Jun would vary more across cell types. 
Because the Jaccard similarty of Jun ChIP-seq of the two given cell types (HepG2 and HeLaS3) is much lower than that of CTCF ChIP-seq.** 

4. Based on these results, do the genomic locations found by ChIP-seq experiments depend more on the cell-type, or on the transcription factor (tf) being assayed? (20 pts)
**Based on the results, the genomic regions found by ChIP-seq experiments depends more on the TF. 
Because the interval sets have highest Jaccard similarty (other than themselves) with interval sets resulted from ChIP-seq of the same TF and different cell type, 
followed by same cell type and different TF. Interval sets of different TF and different cell type have the lowest Jaccard similarty score. 
(except helas3_ctcf: helas3_ctcf vs hepg2_jun has higher Jaccard similarty score than helas3_ctcf vs helas3_jun)**
