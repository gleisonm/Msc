# Resume

Com o controle de qualidade geramos arquivos “limpos” que utilizaremos para a imputação.
Até agora usamos SNPs que estavam presentes apenas no MicroArray, **a imputação é um método estatístico para incrementar a imputação para predizer genotypes que podem nao estar presentes nas amostras**

Essas variantes estão densamente presentes em um painel de referência de phased haplotypes.

O painel de referência contem variantes mais casuais.

Why imputate?
- Increase power
- Fine-mapping
- Meta-analysis

## How it works?
(Talvez fazer um mermaid disso, ou colar as imagens dos slides)
1. Phasing
	1. Set Haplotypes: Prediz estatisticamente quais são os possíveis haplótipos dos pais.
2. **Reference painels**
	1. Haplotypes from family tree: Alinha os haplótipos preditos com o painel referência
3. Imputation
	1. Imput Haplotypes
4. Imputed Genotypes

Funciona como um alinhamento colocando entre os nucleotídeos informações faltantes que estão presentes nos paineis referência. 
## Tools
Shapeit/imput2
Beagle
MACH/Minimac
## Reference panels
HapMap2
HapMap3
1000 Genomes
**Haplotype Ref consortium**
**TOPMed**

## Types of Imputated genotypes
Isso é uma predição, e para avaliar isso, temos algumas métricas:
- Best-guessed genotype (GT) with maximum posterior probability (most probable genotype)
	- Geralmente não é utilizado, porque ele não conta a incerteza da imputação (**só é utilizado quando não há imputação do dado**).
- Genotype probabilities (GPs) of all three genotypes (ref/ref, ref/alt, alt/alt)
- Genotype dosage (DS) **dosage of alternative alelle**: probability of alternative homozygote genotype. Value can be between 0-2

Beagle output file (VCF format:)
(colar aqui o exemplo)

## Imputation Quality
BEAGLE: r²
IMPUTE2: info score
…
- Value of 1: no uncertain in the imputed genotypes
- Value of 0: complete uncertain about genotypes
- Above 0.4
## Imputation QC
- Pre-imputation: Exclude MAF < 1% variants
- Post-imputation: Information measures between 0-1

Wich
Array sufficiente quality
What haplotypes of reference panels
Which association 

## Summary
….