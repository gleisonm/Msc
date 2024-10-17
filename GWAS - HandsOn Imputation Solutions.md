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

The files specified after the -data option indicate the imputed genotype and sample files, whilst the -o option indicates the name of the SNPTEST output file. The -frequentist option specifies the association model to be fitted: choosing “1” indicates an additive model. The option -method specifies the approach used to fit this model: using “expected” fits a regression model with the genotype probabilities. The option -pheno specifies the name of the column in the sample file containing the phenotype to be tested for association. Finally, the option -cov_names indicates the names of the columns in the sample file to be adjusted for in the association analysis (separated by a single space).