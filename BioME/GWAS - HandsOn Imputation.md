# Population Structure in GWAS

Nesta prática iremos avaliar a evidencia da estrutura populacional em um dataset que já passou pelo controle de qualidade.

Usaremos data do **1000 Genomes Project**. Indivíduos não-Europeus serão removidos.

## Prepare 1000 Genomes data
OBS: **Arquivo pesado, pode passar para o próximo passo já com o arquivo convertido**

Iremos baixar o arquivo que contem informação de 629 indivíduos:
```bash
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
```

Converteremos o .vcf para o formato PLINK:
```bash
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
```

Note que o arquivo **ALL.2of4intersection.20100804.genotypes.bim** contem SNPs sem um rs_ID e estão indicados com “.”

Essas informações faltantes não serão um problema para este tutorial. Entretanto, é uma boa prática definir identificadores únicos para esses SNPs
```bash
plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs
```

O próximo passo será fazer um controle de qualidade.

Removendo variantes baseadas no missing genotype data:
```bash
plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS
```

Removendo indivíduos baseado no missing genotype data:
```bash
plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2
```

Removendo variantes baseado no missing genotype data:
```bash
plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3
```

Removendo indivíduos baseado no missing genotype data:
```bash
plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4
```

Removendo variantes baseado no MAF:
```bash
plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5
```

Extraindo do 1000 Genome data, variantes presentes no nosso HapMap clean dataset
```bash
awk '{print$2}' HapMap_3_r3_clean.bim > HapMap_SNPs.txt

plink --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6
```

Faremos o mesmo para variantes presentes no 1000 Genome dataset, extraindo elas do nosso HapMap clean dataset
```bash
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt

plink --bfile HapMap_3_r3_clean --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out HapMap_MDS
```

**Os datasets agora contem o mesmo número de variantes. Esses datasets devem ser da mesma build.**

Criaremos agora o arquivo que contêm um SNP_ID e a posição física por linha
```bash
awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
```

Agora iremos usar essa informação para mudar a build do 1k Genome data:
```bash
plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
```

## Merge the HapMap and 1k Genome datasets
Para fazer esse merge devemos ter certeza que esses arquivos são mergeable:
1. Tenha certeza que **o genoma referência é o mesmo nos dois datasets**
2. Resolva as strand issues
3. Remova os SNPs depois dos dois últimos passos que diferem entre os dois datasets

Os próximos serão para comparar os dois datasets e ter certeza que eles se correspondem:

1. Set o genoma referência:
```bash
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt

plink --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj
```

O 1kG_MDS7 e o HapMap-adj tem o mesmo genoma referencia para todos os SNPs. Este comando ira gerar alguns warnings para impossibles A1 allele assignment

2. Resolver Strand Issues
```bash
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp

awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp

sort 1kGMDS7_tmp HapMap-adj_tmp | uniq -u > all_differences.txt
```

2086 diferenças entre os arquivos, alguns desses pode ser por causa de strand issue.

Iremos dar um flip nos SNPs para resolver esses strands issues. Iremos printar o SNP_ID e remover duplicatas.

[**--flip**] swaps A↔T and C↔G
```bash
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
  
#Flip the 1043 non-corresponding SNPs:

plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap 
```

Vamos checar quais SNPs continuam problemáticos mesmo após o flip.
```bash
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp

sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u > uncorresponding_SNPs.txt
```

Este arquivo mostra que há 136 diferenças entre os arquivos

3. Remover SNPs problemáticos dos dois arquivos
```bash
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
```

O comando acima gera uma lista dos 68 SNPs que são discordantes entre os dois arquivos, por isso iremos remove-los:
```bash
plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS2

plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8
```

Agora iremos juntar os dois arquivos:
```bash
plink --bfile HapMap_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2
```

**Note, we are fully aware of the sample overlap between the HapMap and 1000 Genomes datasets. However, for the purpose of this tutorial this is not important.**

## Peforming multi-dimensional scaling

Multi-dimensional scaling no PLINK é feito em 3 passos:
1. **Itentificar o subset de SNPs comuns independentes** (minor allele frequency at least 5%) que não estão lincados com outros
2. **Calcular o relatedeness entre cada par de indivíuos** usando o set de SNPs independentes que são a medida de similaridade
3. **Performar multidimensional scaling** usando a matriz de correlação

O primeiro passo é identico ao utilizado na prática 1 [[GWAS - HandsOn QC]]

```bash
plink --bfile MDS_merge2 --extract plink.prune.in --Z-genome --out MDS_merge2

plink --bfile MDS_merge2 --read-genome MDS_merge2.genome.gz --cluster --mds-plot 10 --out MDS_merge2
```

Note que o parâmetro [–Z-genome] cria um arquivo compactado do genoma.
O parâmetro [–mds-plot] especifica quantos eigenvectors da multi-dimensional scaling to summarise in the output file.

Os eigenvectors são o output para o arquivo plink.mds. Cada linha é um indivíduo contendo informações dos identificadores no .fam file, junto com os 10 eigenvectors (C1 e C2)

Para gerar o MDS plot, primeiro necessitamos baixar o arquivo com a infomação populacional do 1k Genome dataset.
```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel
```

Agora podemos converter o código das populações em superpopulation codes (i.e., AFR, AMR, ASN and EUR).
```bash
awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt

sed 's/JPT/ASN/g' race_1kG.txt > race_1kG2.txt

sed 's/ASW/AFR/g' race_1kG2.txt > race_1kG3.txt

sed 's/CEU/EUR/g' race_1kG3.txt > race_1kG4.txt

sed 's/CHB/ASN/g' race_1kG4.txt > race_1kG5.txt

sed 's/CHD/ASN/g' race_1kG5.txt > race_1kG6.txt

sed 's/YRI/AFR/g' race_1kG6.txt > race_1kG7.txt

sed 's/LWK/AFR/g' race_1kG7.txt > race_1kG8.txt

sed 's/TSI/EUR/g' race_1kG8.txt > race_1kG9.txt

sed 's/MXL/AMR/g' race_1kG9.txt > race_1kG10.txt

sed 's/GBR/EUR/g' race_1kG10.txt > race_1kG11.txt

sed 's/FIN/EUR/g' race_1kG11.txt > race_1kG12.txt

sed 's/CHS/ASN/g' race_1kG12.txt > race_1kG13.txt

sed 's/PUR/AMR/g' race_1kG13.txt > race_1kG14.txt
```

Crie um arquivo de raça do seu proprio dataset:
```bash
awk '{print$1,$2,"OWN"}' HapMap_MDS.fam > racefile_own.txt

# Concatenate race files:
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt
```

Gere o plot de estratificação:
```R
data<- read.table(file="MDS_merge2.mds",header=TRUE)
race<- read.table(file="racefile.txt",header=TRUE)
datafile<- merge(data,race,by=c("IID","FID"))
head(datafile)

pdf("MDS.pdf",width=7,height=7)
for (i in 1:nrow(datafile))
{
if (datafile[i,14]=="EUR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")}
par(new=T)
if (datafile[i,14]=="ASN") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="red")}
par(new=T)
if (datafile[i,14]=="AMR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col=470)}
par(new=T)
if (datafile[i,14]=="AFR") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="blue")}
par(new=T)
if (datafile[i,14]=="OWN") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")}
par(new=T)
}

abline(v=-0.035,lty=3)
abline(h=0.035,lty=3)
legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red",470,"blue","black"),bty="o",cex=1)
```

![[MDS.pdf]]

Este plot demonstra onde seu aquivo recai sobre o grupo europeu do 1k Genome data. Entretanto não temos outliers de raça.

For educational purposes however, we give scripts below to filter out population stratification outliers. Please execute the script below in order to generate the appropriate files for the next tutorial.

## Exclude Ethnic Outliers
Selecione os indivíduos no HapMap data abaixo do threshold. 