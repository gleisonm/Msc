## Tópicos Principais
- Banco de Dados
- Microarray

## Conceitos-Chave
- 

## Resumo
### Illumina Microarrays
- Table com os diferentes sinais de luz para cada SNP

### Genome Studio
- Genotype calling algorithms
	- Analyse raw intensity data
	- Estimate prob ….

### Cluster Plot
- Observar no clusterplot a separação de cada SNP  (devem estar bem separados)
### QC
#### Per-individual QC
- Gender checks
	- Plating errors
	- Male 1 x ChrX heterozygous
	- **Male DNA samples marked as female in input files will have a higher tan expected hemozigosity rate **
	- **Female samples marked as male will haver a lower than expected heterozygosity rate**
	- Men F of homozygosity > 0.8
	- women F homozigosity < 0.2
	- PLINK gender test
- call rate
	- Individuals with more than 1-5% missing genotypes are removed
	- Sample call rate of 95% is usually used
- Sample heterozygosity
	- missing genes vs heterozygosity
	- Genotipe failure: Genotypes are classified as missing ….
- Identity-by-state (IBS) - PLINK
	- identical samples will share IBS near to 100% (allowing for geotyping errors)
	- Related individuals will share higher IBS than unrelated individuals
	- Identical by descent (IBD) two or more individuals if they …
		- IBD = 1 duplicated or monozygotic twins
		- 0.5 for the first degree ….
	- PI_HAT: high = good sample
		- sample contamination
		- pedegree errors
		- unknown familial relationship
- PCA and MDS
	- Population substructure
	- ….

##### SNP QC (Review)
- SNP call rate
	- Better = call rate > 95%
	- **difference between missing and present SNPs**
- HWE
	- Indicative of:
		- genotype calling error
		- **only for controll**
- MAF
	- Remove individuals with minor allele frequency (?)
- Mendelian errors
	- Identify impossible combinations of alleles due to family tree
- Genome builds, ref genome and variants alignment
	- Position of SNPs change
	- some SNPs change strands between builds
	- Older arrays were generated using previews genome builds (**lift over needes**)
	- lift-over - Change builds

##### Phenotype QC
- Data checklist
	- biologically impossible values
	- missing values
	- Impossible phenotypes and extreme phenotypes (like 500kg)
	- ….

##### Summary
Participants need to provide:
- Ethnic background
- Age distribution
- Gender ratio
- Disease prevalence

#### Post-association QC
1. Allele frequency diferences can indicate….


#class #gwas #microarray