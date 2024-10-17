In this practical, we will combine association summary statistics from the imputed 1Mb region considered in Practical 3 from two studies to increase the strength of the signal using fixed- and random-effects meta-analysis in GWAMA. You will combine the results of your analyses from Practical 4 with results from an additional study.

More details on meta-analysis with GWAMA can be found at the following website:

[https://genomics.ut.ee/en/tools](https://genomics.ut.ee/en/tools)

# Convert SNPTEST output to GWAMA input file format

GWAMA requires one input file per study. Each row of the file provides association summary statistics for a SNP, arranged in columns. For each SNP, the following columns are required: SNP ID, effect allele, other allele, effect size and standard error (or odds ratio and 95% confidence limits). The default column headers for these statistics are: MARKERNAME, EA, NEA, BETA and SE.

You can extract the relevant summary statistics from the SNPTEST output generated in Practical 4 by typing the following commands:

```bash
sed '/^#/d' STUDY.snptest > STUDY.snptest1

awk 'NR==1; NR > 1{if($9>0.4){print $0;}}' STUDY.snptest1 > STUDY.snptest2

perl SNPTEST2_2_GWAMA.pl STUDY.snptest2 gwama.STUDY.stage1.txt SE
```

Change in gwama.STUDY.stage1.txt file MARKER column into MARKERNAME using nano editor.

Above, we have used an imputation info of 0.4, as a cut off for selecting “well imputed” variants.

On completion, the program will create the GWAMA input file gwama.STUDY.stage1.txt.

# Additional files required by GWAMA

We have created a GWAMA input file for a second study, containing association summary statistics for the same 1Mb region, and filtered for imputation info of 0.4. The file is named gwama.STUDY.stage2.txt, and has similar format as the file you have created in this practical.

GWAMA requires an additional file, by default gwama.in, which lists the files to be combined in meta-analysis. You should use a text editor to create a file called **gwama.in**, and type the filename of each study to be included in the meta-analysis, in this case gwama.STUDY.stage1.txt and gwama.STUDY.stage2.txt, one on each row.

# Fixed-effects meta-analysis in GWAMA

GWAMA implements fixed-effects meta-analysis of association summary statistics for all studies listed in the gwama.in file using an inverse-variance weighting scheme. You can perform the meta-analysis in GWAMA using the following command:

  

./GWAMA -qt -o gwama.fixed

  

The -qt flag tells GWAMA that the association summary statistics are stored as effect size and standard error (and no odds ratio and 95% confidence limits). The -o flag determines the name of the GWAMA output file.

  

Take a look at the first few lines of the output file using the command:

  

head gwama.fixed.out

  

The file contains one row per SNP, and provides the following information: SNP ID, reference (effect) allele and other allele, effect allele frequency (-9 for missing here as frequencies not provided in input files), meta-analysis effect size and standard error, 95% confidence limits, Z-score, p-value and –log10 p-value, Cochran’s Q statistic and p-value for heterogeneity, I2, number of studies and samples for which the SNP is reported, and the direction of effects in each contributing study (aligned to reference allele, +, -, 0, or ? if missing).

  

To view the SNPs with the strongest signal of association after meta-analysis, you can use the command:

  

sort -g -k10 gwama.fixed.out > gwama.fixed.sort

  

head gwama.fixed.sort

  

Which SNP shows the strongest signal of association? What is the rsID and p-value? Is there any evidence of heterogeneity in allelic effects of the SNP between the two studies?