## Files
1. [HeapMap.bim](https://www.cog-genomics.org/plink/1.9/formats#bim) 
	- Cada linha é um SNP, com informação sobre o cromossomo e posição
2. [HeapMap.fam](https://www.cog-genomics.org/plink/1.9/formats#fam).
	- Cada linha é um indivíduo, contendo o id, parental_id, sex e status da doença
1. [HeapMap.bed](https://www.cog-genomics.org/plink/1.9/formats#bed) 
	- É o link entre esse dois arquivos.

Estes arquivos estão conectados entre si, **caso necessite remover alguma informação, é necessário remover de todos os arquivos**

Please note: There are other formats on input files, such as VCF, oxford format (gen/bgen)…

![[Pasted image 20241014132241.png]]

Em alguns casos podem haver informações faltantes nos arquivos, para contornar isso é possível acrescentar os parâmetros seguintes:
```bash
--no-fid
--no-parents
--no-sex
--no-pheno
```

Genótipos faltantes geralmente são representados por **0** no arquivo **.ped** ou arquivos similares. É possível mudar isso para qualquer outro caracter utilizando o parâmetro [–missing-genotype], **entretanto 1, 2, 3, 4, A, C, G, T** não é permitido.

Em arquivos .bim “.” também é considerado como valor NA (**isso só foi permitido após 16 jan 2023**), esse valor pode ser alterado com [–missing-genotype2].

Com o PLINK é possível gerar um summary das informações rodando o comando:
```bash
plink --bfile HapMap_3_r3_1 --make-bed
```

## Sample QC

### Primeiro Passo - Low Quality Samples

O primeiro passo do QC é identificar amostras com qualidade baixa, as principais métricas utilizadas são:
- cell rate
	- Low cell rates indicam de baixa qualidade do DNA (Untrustable)
- heterozygosity
	- Baixos níveis de heterozygositu são indicativos de inbreeding (?)

You can use PLINK to generate summaries of these metrics using the following command:
```bash
plink --bfile HapMap_3_r3_1 --missing --het --autosome --out plink.sampleqc
```

Este comando **cria dois arquivos relevantes para o QC**:
1. plink.sampleqc.imiss
	- Entrega um resumo dos números de missing genotypes para cada indivíduo.
	- Cada linha é um indivíduo contendo informações de número (N_MISS) e proporção (F_MISS) de missing genotypes
1. plink.sampleqc.het
	- Entrega um resumo do número de hemozygous genotype calls para cada indivíduo.
	- Cada linha é um indivíduo contendo informações do número de homozygous genotypes calls (HOM), e o total de number of genotype calls (NM)

Não há uma pré-determinação de tresholds para o controle de qualidade, tudo irá depender da fonte do DNA, o array utilizado, a distribuição da frequencia alélica e outros fatores.

**Sample QC foca em encontrar outliers que podem facilmente serem identificados por gráficos das métricas vs SNPs**

Scatterplot of sample cell rate against heterozygosity
```R
indmiss<-read.table(file="plink.sampleqc.imiss", header=TRUE)
snpmiss<-read.table(file="plink.sampleqc.lmiss", header=TRUE)
het<-read.table(file="plink.sampleqc.het", header=TRUE)
d3<-merge(indmiss,het, by="IID")
d3$het<-(d3$N.NM.-d3$O.HOM.)/d3$N.NM.

# read data into R 

pdf("histimiss.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness") #selects column 6, names header of file

pdf("histlmiss.pdf") 
hist(snpmiss[,5],main="Histogram SNP missingness")  

pdf("imiss_het_plot.pdf") 
plot(d3$het, 1-d3$F_MISS, xlab="Heterozygosity", ylab="Callrate")
```

![[histimiss.pdf]]

Cada ponto corresponde à uma amostra, de acordo com seu call rate (eixo-y) e sua heterozygosity (eixo-x).

The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean:
```R
het <- read.table("plink.sampleqc.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)
```

![[histlmiss.pdf]]

De forma alternativa, podemos gerar uma lista de indivíduos que não atingem os critérios de qualidade:
```R
callrate <- as.double(0.9)
heter_low <- as.double(0.305)
heter_upp <- as.double(0.4)
d1<-read.table(file="plink.sampleqc.imiss", header=TRUE)
d2<-read.table(file="plink.sampleqc.het", header=TRUE)
d3<-merge(d1,d2, by="IID")
d3$het<-(d3$N.NM.-d3$O.HOM.)/d3$N.NM.
write.table(d3[(!(1-d3$F_MISS > callrate & d3$het < heter_upp & d3$het > heter_low)),2:1], file= "", quote=F, col.names = F, row.names = F)
```

Com esse comando, criamos uma lista de indivíduos a serem excluidos com base no call rate < 0.9 e heterozygosity < 0.305 ou > 0.4

Para remover esses indivíduos da análise iremos gerar um novo dataset chamado **HapMap_3_r3_2** com o seguinte comando:

```bash
plink --bfile HapMap_3_r3_1 --remove fail-imiss-het-qc.txt --make-bed --out HapMap_3_r3_2
```

De forma alternativa podemos utilizar o parâmetro [–min] para remover indivíduos com missing call rates excedendo o valor definido, mais info [aqui](https://www.cog-genomics.org/plink/1.9/filter#missing)

### Segundo Passo - Sex Disparity

O segundo passo será checar a discrepância pelo gênero. Para isso identificaremos **plating erros e mistura de amostras** utilizando as infomações do genótipo pelo cromossomo X.

Premissas:
- Homens só possuem uma cópia do cromossomo X
	- Por isso não podem ser heterozigoto para nenhum marcador que não está na região pseudoautossomal(?) do chrY.
	- Quando um algorítimo de genotyping calling detecta a male heterozygote para um marcador do chrX ele informa que o genotypo é missing.
- Mulheres que são **erroneamente** identificadas como Homens  terão grandes quantidades de missing data porque todo o chrX será setado como missing.
	- Terão menos do que esperado homozygosity rate
- Homens que são **erroneamente** identificados como mulheres terão mais do que esperado homozygositi rate

Nem todo algorítimo fará essa conversão. **Typically, one expects male samples to have a homozygosity rate of 1 (although, because of genotyping error, this can vary) and females to have a homozygosity rate of <0.2**.

A melhor forma de identificar discrepâncias relacionadas ao gênero é calcular o média da homozygosity rate (Xchr inbreeding) across X-chr markers for each individual in the study.

**Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK**

Para fazer QC de acordo com o genero rodaremos o seguinte comando:
```bash
plink --bfile HapMap_3_r3_2 --check-sex
```

Podemos visualizar o resultado do sex-check usando o seguinte script:

```R
gender <- read.table("plink.sexcheck", header=T,as.is=T)

pdf("Gender_check.pdf")
hist(gender[,6],main="Gender", xlab="F")
dev.off()

pdf("Men_check.pdf")
male=subset(gender, gender$PEDSEX==1)
hist(male[,6],main="Men",xlab="F")
dev.off()

pdf("Women_check.pdf")
female=subset(gender, gender$PEDSEX==2)
hist(female[,6],main="Women",xlab="F")
```

![[Gender_check.pdf]] Gender

![[Men_check.pdf]] Man


![[Women_check.pdf]] Woman

Isso indica que existe uma mulher com discrepância de gênero, F value = 0.99

Iremos criar uma lista de indivíduos que estão com o gênero descordando:
```bash
grep "PROBLEM" plink.sexcheck > fail-sexcheck-qc.txt
```

A coluna 1 é a family ID e a coluna 2 é o indivíduo ID, a coluna 3 ajusta o sexo e a coluna 4 denota o sexo de acordo com genotype data.

**When the homozygosity rate is more than 0.2 but less than 0.8, the genotype data are inconclusive regarding the sex of an individual and these are marked in column 4 with a 0.**

Gerando lista de exclusão:
```bash
grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}'> fail-sexcheck-qc.txt
```

Para remover esses indivíduos da análise, e gerar um novo dataset chamado HapMap_3_r3_3, usaremos o seguinte comando:
```bash
plink --bfile HapMap_3_r3_2 --remove fail-sexcheck-qc.txt --make-bed --out HapMap_3_r3_3
```

### Terceiro Passo - SNP QC

Agora iremos identificar SNPs com baixa qualidade.
As métricas mais utilizadas são:
 - call rate
	 - Baixos call rates indicam baixa qualidade
 - Deviation from HWE
	 - Altos desvios indicam baixa qualidade


