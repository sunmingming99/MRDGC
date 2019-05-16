# MRDGC
MRDGC: a parallel approach for the identification of master regulators based on the differentially expressed genes and the regulatory capacity of regulators
====================================================


Data preparation
----

You need to prepare the expression profiles of genes, transcription factors and microRNA.
    For example, you can download all the expression profiles from The Cancer Genome Atlas (TCGA) (https://cancergenome.nih.gov/) data portal. And then filter expression profiles of transcription factors from TCGA and The Human Transcription Factors database (http://humantfs.ccbr.utoronto.ca/download.php).
    The folder named "Experimental Data" lists the data we used in the article as described in section "Data Collection and Pre-process".
    "GENE_EXP.txt" is the expression profile of differentially expressed genes which ranked descending according |logFC| after analysing the differential expression of genes. "TF_EXP.txt" and "miRNA_EXP.txt" are the expression profiles of transcription factors and miRNAs, espectively.


Parameters
----

threshold theta=0.01
TOP K=100



Usage
----

Parrllel.m is the parallel implementation of the MRDGC in Matlab.
    run Parrllel.m


    "regulators_pvalue.csv" records the p-value of all regulators;
    res_candite_totalgenes.csv" records the candidate master regulators;
    "master_regultors.csv" records the final master regulators which we select TOP K from candidate master regulators.


Sequence.m is the serial implementation of MRDGC in Matlab which we used to compare with parrllel algorithm.
    run Sequence.m, the result is also recorded in "regulators_pvalue.csv".


    "regulators_pvalue.csv" records the p-value of all regulators;
    res_candite_totalgenes.csv" records the candidate master regulators;
    "master_regultors.csv" records the final master regulators which we select TOP K from candidate master regulators.


Result
----

The result is listed in "master_regultors.csv" when we set the threshold of p-value (pvalue<0.01) and seclte TOP K (K=100) regulators as master regulators. 




Mingming Sun, mm_sun@hnu.edu.cn
2019.5.3
