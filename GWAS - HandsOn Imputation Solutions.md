Nessa prática iremos permorf a imputation of a 1Mb region flanking de lead SNP from association analysis [[GWAS - HandsOn Association Analysis]] na finalidade de aumentar o sinal e fine-map a variante.

Iremos utilizar a pre-phased scaffold gerado de um dataset limpo. Imputaremos up to the março de 2012 [todos ancestrais] reference painel do 1k Genome Project.

Utilizaremos o SNPTEST e EPACTS para testar a análise de associação para todos os SNPs desta região a partir de um modelo aditivo, com o ajuste de eigenvectors de uma escala multi-dimensional para account for population structure.

More details on **pre-phasing with EAGLE** and **imputation with BEAGLE** can be found at the following websites:

- [Eagle](https://data.broadinstitute.org/alkesgroup/Eagle/)
- [Beagle](http://faculty.washington.edu/browning/beagle/beagle.html)

**BEAGLE website also contain links to downloadable formatted reference files required for imputation up to reference panels from the 1000 Genomes Project.**

# Generating a clean genotype scaffold and pre-phasing with EAGLE

Passos para Imputation:
1. Gerar um scaffold limpo de genotype data para cada cromossomo baseado no processo de controle de qualidade [[GWAS - HandsOn QC]]
2. Pre-Phase the scaffold para cada cromossomo separadamente utilizando EAGLE e reference gentic maps refletinfo o rate de recombinação (avaliable from BEAGLE)
**Esse processo é computacionalmente custoso, entretanto faremos isso apenas para uma  região de 1Mb**

## Cleaning scaffold for chromossome 1
- Removeremos SNPs com MAF < 1% (They are more prone to genotyping errors)
- Using cleaned file from [[GWAS - HandsOn Imputation]]:
```bash
plink --noweb --bfile STUDY.clean --maf 0.01 --chr 1 --recode vcf-iid --from-kb 75629 --to-kb 76629 --out STUDY.scaffold.chr1
```

Agora utilizaremos o EAGLE par pre-phase the scaffold para o cromossomo1:
```bash
/path/Eagle_v2.4.1/eagle --vcf STUDY.scaffold.chr1.vcf --geneticMapFile genetic_map_chr1_b37.txt --numThreads=2 --Kpbwt=20000 --outPrefix STUDY1_prephased
```

Neste comando, o parâmetro [–geneticMapFile] especifica o arquivo contendo rates de recombinação across the cromossome, e pode ser baixado com o eagle.

The genotype file in vcf format contains several header rows, each starting with “#”. They contain information about vcf format version (v.4.2 in our case). The last header line contains following columns:

- CHROM – chromosome

- POS – position

- ID – marker ID

- REF – reference allele (please note that this doesn’t have to be the same as human genome reference as noted in header)

- ALT – alternative allele

- QUAL – information about marker quality (no information is marker as “.”)

- FILTER – information if the marker passed the quality filters

- INFO – marker annotation (currently only containing “PR”, which means provisional reference allele)

- FORMAT - It defines, how the genotype data is coded in given file. In current file, the genotypes are defined as: 0/0 is reference allele homozygote; 0/1 is heterozygote, and 1/1 is homozygote of alternative allele.

Following columns contain sample names and their genotype data in following rows.

VCF format is a standard format for storing whole genome and exome sequencing data. You can find more information about VCFv4.2 from:

[VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

Most commonly used tools for manipulating VCF formatted files are VCFTOOLS and BCFTOOLS software tools. VCF is one of the standard formats, used by many genetic data analysis tools.

# Imputation of pre-phased scaffold
A imputação pode ser dividida entre non-overlapings chunks (geralmente 1Mb) através dos cromossomos, incluindo um buffer por volta de 500kb em ambos os lados para evitar edge effects.

```bash
beagle gt=STUDY1_prephased.vcf.gz ref=1000G_ref_chr1_75629223-76629223.vcf.gz map=plink.chr1.GRCh37.map out=STUDY_imputed nthreads=4 ne=20000 gp=true
```

In the above command, the gt, ref, and map options specify the phased genotype scaffold data, the reference panel haplotype file for 75629223-76629223 region of the chromosome 1, and the genetic map file for chromosome 1 (as used by EAGLE but in slightly different format), respectively. Both reference data as well as genetic map files can be downloaded from the BEAGLE website. The _nthreads_ option defines the number of computing cores used for imputation. Other options specify the output file (_out_), that we want to get genotype probabilities in the imputation output file (_gp_), the effective sample size (_ne_, 20000 is a default recommended by the BEAGLE developers).

The file format is again similar to those of the previous genotype data files. There are a few additions through: INFO column now includes information about imputation quality (DR2), allele frequency of imputed variants (AF), and IMP value indicating that the variant is imputed. Genotype data now has three different fields separated by “:”. These three fields are:  
GT: most probable genotype

DS: dosage of alternative allele, which is calculated as probability of heterozygous genotype + 2 * probability of alternative homozygote genotype. Value can be between 0-2.

GP: probabilities of all three genotypes.

For example, genotype value “0|1:0.57:0.43,0.57,0”, means that the most probable genotype in this position for observed sample is heterozygote (marked as 0|1), the dosage of alternative allele is 0.57 (=0.57+2*0), and the probabilities of the genotypes are: REF/REF 43%; REF/ALT 57%; ALT/ALT 0%. It means that even though the probability of heterozygote is the highest, it can in fact be reference homozygote and the imputation certainty is not very good in this particular case.

Let’s check another case with a genotype value “0|0:0:1,0,0”. In this case, the most probable genotype is homozygote reference (0|0), the dosage of alternative allele is zero and genotype probabilities are: REF/REF 100%; REF/ALT 0%; ALT/ALT 0%. Therefore, the imputation certainty for this genotype is good – it has a very high probability of being reference allele homozygote.

Based on the certainty of genotypes and the variant minor allele frequency, the imputation quality of a variant is estimated (DR2), where the value above 0.4 indicates a well imputed variant. Usually, this threshold is used for variant filtering in large GWAS consortia.

# Testing for association using SNPTEST
Before running SNPTEST, we have to convert the files from the VCF format into IMPUTE format. Even though SNPTEST should work with VCF format files, the current version gives an error, while trying to run it with BEAGLE imputed files. Therefore, first we will remove multiallelic variants using bcftools and then we will use PLINK2 to convert our genotype data:  

```bash
bcftools index STUDY_imputed.vcf.gz

bcftools view --max-alleles 2 --exclude-types indels STUDY_imputed.vcf.gz -o STUDY_imputed_nomultiallelic.vcf

plink2 --vcf STUDY_imputed_nomultiallelic.vcf dosage=GP-force --recode oxford --out STUDY.impute --double-id
```

We will next use SNPTEST to test for association of SNPs in the 1Mb imputed region under an additive model, adjusting for eigenvectors from multi-dimensional scaling to account for population structure. To perform this analysis, SNPTEST requires two files: 
1. an imputed genotype file, generated by IMPUTE;
2. sample file, created based on the phenotype file.
	- The sample file should contain one row per individual, with two headers as defined above. 
	- The sample file should contain information on the phenotype to be analysed, together with any relevant covariates that will be adjusted for in the analysis. 
	- In the sample file, a binary phenotype should be coded as “1” for cases and “0” for controls (in contrast to PLINK, which codes cases and controls as “2” and “1”, respectively). 

To save time for you in this practical, we have already created the *.sample file. We have named the sample file COHORT.impute.sample

```bash
#You can now download SNPTEST:
wget [http://www.well.ox.ac.uk/~gav/resources/snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz](http://www.well.ox.ac.uk/~gav/resources/snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz)

tar -zxvf snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz
```

You can test for association with the phenotype (pheno) under an additive model, adjusting for two eigenvectors (C1 and C2), using the following SNPTEST command (type on one continuous line):

```bash
snptest -data STUDY.impute.gen COHORT.impute.sample -o STUDY.snptest -frequentist 1 -method expected -pheno pheno -cov_names C1 C2
```

The files specified after the -data option indicate the imputed genotype and sample files, whilst the -o option indicates the name of the SNPTEST output file. The [-frequentist] option specifies the association model to be fitted: choosing “1” indicates an additive model. The option -method specifies the approach used to fit this model: using “expected” fits a regression model with the genotype probabilities. The option [-pheno] specifies the name of the column in the sample file containing the phenotype to be tested for association. Finally, the option [-cov_names] indicates the names of the columns in the sample file to be adjusted for in the association analysis (separated by a single space).

The command will take some time to run, so please be patient. Output will be printed to the screen, including summaries of the phenotype and covariates to be included in the analysis. The program analyses data in “chunks” of 100 SNPs, and will report progress to the screen. When SNPTEST finishes running, you should see the message “finito”.

The file is very detailed and contains a lot of information! After some initial information confirming the association test used, the file contains one row per SNP, with columns providing different summary statistics for the quality of imputation and the association with the phenotype. The most important columns in the output are:

- rsid: identifier of the SNP (in the 1000 Genomes reference panel)
- position: position of the SNP (basepairs)
- alleleA and alleleB: the pair of alleles at the SNP
- info: the “info score” metric of imputation quality
- frequentist_add_pvalue: p-value for association from the additive model
- frequentist_add_beta_1: the log-odds ratio for allele B relative to allele A
- frequentist_add_se_1: the standard error of the log-odds ratio
- comment: information regarding problems with fitting the regression model

For many of the SNPs, you will see a comment about the model not being fitted, and many of the association summary statistics are “NA”. This occurs typically for rare variants, where the regression model cannot be fitted because the counts of the minor allele are so low in cases and/or controls.  

You can plot a summary of association statistics from the additive model across the region for all SNPs with info score above a pre-specified threshold using the command:
```R
args<-commandArgs(trailingOnly=TRUE)
name <- args[1]
name2 <- args[2]
n <- as.double(args[3])
out <- args[4]

data<-read.table(name, header=T, sep=" ")

library(vcfR)
library(plyr)

data <- read.table(name, header=T, sep=" ")
d <- read.vcfR(name2)

d_imp <- as.data.frame(d@fix[,c("ID", "INFO")], stringsAsFactors=F)

d_imp <- within(d_imp, {
 rsid <- ID
 imputed <- sapply(strsplit(INFO, ";"), "[", 3)
})

data <- join(data, d_imp, by="rsid", type="l")

data$imputed <- ifelse(data$imputed %in% "IMP", 1, 0)


pdf(out)
plot(data$position[data$info>=n],-log10(data$frequentist_add_pvalue[data$info>=n]), pch=20, xlab="Position (Bp)", ylab="-log10(p-value)",col=c("red","black")[data$imputed+1])
points(data$position[data$info>=n & data$imputed==0],-log10(data$frequentist_add_pvalue[data$info>=n & data$imputed==0]), pch=20,col="red")
```


**IMAGEM**

n this plot, each point represents a SNP, plotted according to the physical position on the x-axis and the -log10 p-value for association from the additive model on the y-axis. SNPs are coloured in red if directly typed, and coloured in black if imputed.

Is the strongest signal of association from a directly typed or imputed SNP?

You’ll notice that all the most strongly associated SNPs localise to the same region. Why do you think this is? To identify which SNPs attain genome-wide significance, you can use the command:

```R
#R --vanilla --slave --args STUDY.snptest 0.4 5e-8 < snptest.extract.R | sort -k 8 -g
R --vanilla --slave --args STUDY.snptest 0.4 5e-8 < snptest.extract.R | sort -k 8 -g
```


How many SNPs attain genome-wide significance in this analysis? Which SNP shows the strongest signal for association (the lead SNP)? What is the log-odds ratio and the corresponding standard error for the lead SNP?

# Optional: Testing for association using EPACTS

EPACTS software is a fast tool for association analysis that enables several statistical tests.

As the first step, we will create an index for our imputed genotype file using tabix software:

```bash
tabix-0.2.6/tabix -p vcf STUDY_imputed.vcf.gz
```

This makes possible to very quickly retrieve data from any position of the file and is necessary to run EPACTS analysis.

This is the sample file, we need to create for running association analysis on EPACTS.

You can test for association with the phenotype (pheno) under an additive model, adjusting for two eigenvectors (C1 and C2), using the following EPACTS command (type on one continuous line):

```bash
epacts single --vcf STUDY_imputed.vcf.gz --ped STUDY.ped --chr 1 --pheno pheno --cov C1 --cov C2 --test b.score --out STUDY_imputed --field DS --run 1
```
  

The file specified after the --vcf option indicate the imputed genotype, and --ped phenotype file, whilst the --out option indicates the name of the EPACTS output files. The --test option specifies the association test, which in our case is logistic score test using additive genetic model. The option --chr indicates that we are analysing chromosome 1, --field tells which genotype data field to use from the vcf file (options are GT, DS, and GP in our case), and --run tells how many cores to use for the calculation. If the --run parameter is not set, EPACTS software generates necessary files for the association analysis but doesn’t run tests themselves. The option --pheno specifies the name of the column in the sample file containing the phenotype to be tested for association. Finally, the option --cov indicates the names of the columns in the sample file to be adjusted for in the association analysis (each covariate with additional --cov argument).

EPACTS generates several output files. The most important ones are:

- STUDY_imputed.epacts.gz - Main results file for the association analysis. File is gzip compressed.
- STUDY_imputed.epacts.mh.pdf - Manhattan plot of the results
- STUDY_imputed.epacts.qq.pdf - Quantile-quantile plot of the results
- STUDY_imputed.epacts.top5000 – Top 5000 associations sorted by the p-value.

The most important columns of the output file are:

- CHROM: chromosome
- BEGIN: position
- MARKER_ID: marker name
- NS : Number of phenotyped samples with non-missing genotypes
- AC : Total Non-reference Allele Count
- MAF : Minor allele frequencies
- PVALUE : P-value of single variant test
- AF.CASE : Non-reference allele frequencies for cases
- AF.CTRL : Non-reference allele frequencies for controls
- 
Download Manhattan and QQ plots for inspection.

How many SNPs attain genome-wide significance in this analysis? Which SNP shows the strongest signal for association (the first SNP in top5000 file).
You can check these by typing

```bash
awk ‘{if($9<5e-08) print $0}’ STUDY_imputed.epacts.top5000 | wc -l

head STUDY_imputed.epacts.top5000
```

How does the association p-values match with the ones from SNPTEST analysis?