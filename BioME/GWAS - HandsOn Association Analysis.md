Essa prática tem conceitos básicos para realizar uma análise de associação com de GWAS em uma regressão logística com o PLINK.

# Testing For Association Under an Additive Model

Para dados de **caso-controle** é comum utilizar a logistic regresssion framework. typically **under the assumption of an additive model** in the effect of the minor allele at a SNP on the log-odds ratio.

Esse modelo é implementado no PLINK e pode ser rodado com o seguinte comando:
```bash
plink --bfile STUDY.clean --logistic --ci 0.95 --out plink.additive
```

O parâmetro [–ci 0.95] especifica que o PLINK deve produzir 95% de intervalo de confiança.

O output file tem uma linha por SNP com informações da posição, o alelo afetado (A1, que por padrão é o minor allele), o teste performado (aqui ADD para aditivo), o número de genótipos chamados por SNP (NMISS), o odds ratio (OR), o desvio padrão do log-odds ratio (SE), o lower e upper 95% limite de confiança (L95 e U95), a associação estatística (STAT), e o p-valor (P)

Podemos plotar a associação estatística deste modelo por cromossomo:
```R
args<-commandArgs(trailingOnly=TRUE)

file=args[1]
criteria=args[2]
outfile=args[3]

plotdata = read.table(file,header=T,sep="",stringsAsFactors=F)
addmodel=plotdata[which(plotdata$TEST==criteria),]

logPval = -log(addmodel$P,base=10)

# plot max and min

max.pos = max(addmodel$BP)
min.pos = min(addmodel$BP)
size.pos = max.pos - min.pos
center.pos = min.pos + ( size.pos / 2 )

# draw plot
pdf(outfile)
plot(addmodel$BP, logPval, xlab="Position", ylab="-log10 p-value", type="p",col=ifelse(logPval>7.30103, "red", "blue"), pch=20, xaxt="n")
axis(1, at=c(seq(from=min.pos,to=max.pos,by=(size.pos/10))))

# write to file

dev.off()
```

**IMAGEM**

Neste plot, cada ponto representa um SNP de acordo com seu locus no eixo x e o -log10 p-valor associado do modelo aditivo no eixo y. **SNPs que são GW significativos estão em vermelho**

To identify which SNPs attain a given significance threshold, you can use the command (on one continuous line):

  

R --vanilla --slave --args plink.additive.assoc.logistic ADD _alpha_ < findalpha.R > gw_significant_results_additive.txt

  

In this command, you replace _alpha_ by the p-value below which you wish to identify SNPs. For example, to identify SNPs attaining genome-wide significance, you would use the command (on one continuous line):

  

R --vanilla --slave --args plink.additive.assoc.logistic ADD 5e-8 < findalpha.R > gw_significant_results_additive.txt

  

  

more gw_significant_results_additive.txt

  

How many SNPs attain genome-wide significance in this analysis? Which SNP shows the strongest signal for association (the lead SNP)? What is the odds ratio for the minor allele and the 95% confidence interval for the lead SNP?

  

  

1. **Testing for association under a general genotypic model**
    

  

An alternative approach to assuming an additive effect of the minor allele on the log-odds ratio is to allow **a more general genotypic model**. This model includes an additional parameter that allows for a deviation from additivity (referred to as a dominance effect). This model can be implemented in PLINK using the following command:

  

./plink --bfile STUDY.clean --maf 0.05 --logistic genotypic --ci 0.95 --out plink.genotypic

  

Have a look at the format of the output file using the command:

  

head plink.genotypic.assoc.logistic

  

The output file now has three rows per SNP. The first row, labelled ADD, corresponds to the additive component of the genotypic model. The second row, labelled DOMDEV, corresponds to the dominance component of the genotypic model, and measures the deviation from additivity. The third row, labelled GENO_2DF, provides the test statistic and p-value for the genotypic model that incorporates both additive and dominance effects. Please ask if you are unsure of any of the output.

  

You can plot a summary of association statistics from the genotypic model across the chromosome using the command:

  

R --vanilla --slave --args plink.genotypic.assoc.logistic GENO_2DF regPLOT_genotypic.pdf < regplot.R

  

In this plot, each point represents a SNP, plotted according to the physical position on the x-axis and the -log10 p-value for association from the general genotypic model on the y-axis. Do any SNPs attain genome-wide significance in this analysis?

  

To identify which SNPs attain genome-wide significance under the genotypic model, you can use the command:

  

R --vanilla --slave --args plink.genotypic.assoc.logistic GENO_2DF 5e-8 < findalpha.R > gw_significant_results_genotypic.txt

  

How many SNPs attain genome-wide significance in this analysis? Which SNP shows the strongest signal for association under the genotypic model?

  

Take a look at the association summary statistics for the lead SNP from your previous analysis under the additive model. You can do this using the command:

  

grep -w _SNPID_ plink.genotypic.assoc.logistic

  

In this command, you replace _SNPID_ with the identifier of the lead SNP. Do you see a stronger signal of association with the general genotypic model? Is there any evidence for a deviation from additivity?

  

  

2. **Accounting for confounders as covariates in a logistic regression model**
    

  

One primary advantage of the logistic regression model is that it is straightforward to **take account of potential confounders (such as sex or age) as covariates**.

  

The PLINK family file (STUDY.clean.fam) includes information about sex in the fifth column. We can then easily test for association under an additive model, accounting for sex as a confounder using the command:

  

./plink --bfile STUDY.clean --logistic --sex --ci 0.95 --out plink.additive.sex

  

This command adds an indicator variable for sex as a covariate in the regression model, coded as 1 for males and 0 for females. Have a look at the format of the output file using the command:

  

head plink.additive.sex.assoc.logistic

  

The output file now has two rows per SNP, the first providing the evidence of association of the SNP with case-control status under an additive model (ADD), and the second providing the evidence of association of sex with case-control status (SEX).

  

Take a look at the association summary statistics for the SNP with the strongest signal of association from your original analysis under the additive model. You can do this using the command:

  

grep -w _SNPID_ plink.additive.sex.assoc.logistic

  

As before, in this command, you replace _SNPID_ with the identifier of the lead SNP. Is there any evidence of association of sex with the outcome? Does inclusion of sex have any impact on the additive effect of the SNP?

  

The inclusion of sex as a covariate is a special case because the information is already recorded in the PLINK family file. To take account of additional confounders requires an additional data file containing the values of all covariates for each individual. The file covar.txt contains additional confounders that we might wish to take account of in the analysis. Take a look at the format of this covariate file using the command:

  

head covar.txt

  

The file contains one row per individual, each with four columns. The first two columns give the identifier of the individual as used in the PLINK family file. The next two columns provide the AGE and BMI of each individual. We can then test for association under an additive model, accounting for age as a confounder using the command (type on one continuous line):

  

./plink --bfile STUDY.clean --logistic --covar covar.txt --covar-name AGE --ci 0.95 --out plink.additive.age

  

Have a look at the format of the output file using the command:

  

head plink.additive.age.assoc.logistic

  

The output file now has two rows per SNP, the first providing the evidence of association of the SNP with case-control status under an additive model (ADD), and the second providing the evidence of association of age with case-control status (AGE).

  

Take a look at the association summary statistics for the SNP with the strongest signal of association from your original analysis under the additive model. You can do this using the command:

  

grep -w _SNPID_ plink.additive.age.assoc.logistic

  

As before, in this command, you replace _SNPID_ with the identifier of the lead SNP. Is there any evidence of association of age with the outcome? Does inclusion of age have any impact on the additive effect of the SNP?

  

We can also test for association under an additive model, with adjustment for multiple confounders as covariates. So, for example, to adjust for sex, age and BMI, we can use the command (type on one continuous line):

  

./plink --bfile STUDY.clean --logistic --sex --covar covar.txt --covar-name AGE,BMI --ci 0.95 --out plink.additive.covars

  

Have a look at the format of the output file using the command:

  

head plink.additive.covars.assoc.logistic

  

The output file now has four rows per SNP, the first providing the evidence of association of the SNP with case-control status under an additive model (ADD), the second providing the evidence of association of age with case-control status (AGE), the third providing the evidence of association of BMI with case-control status (BMI), and the fourth providing the evidence of association of sex with case-control status (SEX).

  

Take a look at the association summary statistics for the SNP with the strongest signal of association from your original analysis under the additive model. You can do this using the command:

  

grep -w SNPID plink.additive.covars.assoc.logistic

  

As before, in this command, you replace SNPID with the identifier of the lead SNP. Is there any evidence of association of the covariates with the outcome? Does adjustment for these covariates have any impact on the additive effect of the SNP?

  

As a final note, sometimes also your phenotype is not included in the plink files but is in an external phenotype file (could be in the same file with the covariates and this is usually the easiest option). Analogous to the --covar and--covar-name, you would then need to specify the --pheno and --pheno-name to indicate the file name in which your phenotype is and the column name for your phenotype.

  

3. **Exploring the results using web-based tools**
    

  

**LocusZoom.** Go to [http://locuszoom.sph.umich.edu/locuszoom/](http://locuszoom.sph.umich.edu/locuszoom/) and click the link Single Plot – Your Data – Original Locuszoom. Make a regional plot using P-values from the association analysis from exercise 1 (plink.additive.assoc.logistic). You should now be familiar with your results file, i.e. which columns report the P-value and Marker names, what is used as a separator, what chromosome and which positions in Mb does the results file cover.

  

Now plot the region on Chr 1 from 75 to 77 Mb. Select hg19/1000 Genomes Mar 2012 EUR as the Genome Build.

  

Study your plot to answer the following questions:

What does each dot represent in the figure?

Which is the most statistically significantly associated SNP in the region?

Which gene is the closest to this SNP/in which gene is this SNP located?

What do the different colours tell us?

  

**NCBI dbSNP.** Next, go to [http://www.ncbi.nlm.nih.gov/SNP/](http://www.ncbi.nlm.nih.gov/SNP/) and type in as a search term your most statistically significant SNP rsid.

  

What does the website tell us about this SNP?

For example, what is its MAF? In which gene does it belong to?

  

**UCSC genome browser.** Let’s continue our investigations with the UCSC genome browser: [https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu)

Type in the SNP rsid and hit go.

  

Play with zooming of the region and the options underneath the plot.

  

**1000 Genomes.** The 1000 Genomes Project has also several genome browsers: [http://www.1000genomes.org/1000-genomes-browsers](http://www.1000genomes.org/1000-genomes-browsers)

  

Let’s see what the [**Ensembl GRCh37**](http://grch37.ensembl.org/Homo_sapiens/Info/Index) says about our top variant. Explore the variation using the different databases available, e.g. genomic context and phenotype data.

  

**GWAS Catalog.** Finally, let’s see what GWAS catalogue can tell us about our top variant/gene. Go to [https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/) and again, type in the rsid and/or gene.

  

Which traits/diseases has the SNP/gene been previously associated with?

  

**LDlink.** Go to https://ldlink.nih.gov/?tab=home

  

Imagine we would like to replicate our result in independent samples but our collaborators would not have our top SNP genotyped. Which would be the best proxies for our SNP? Use the LDproxy sheet for these calculations. Remember that it is important to select the correct population.

  

Sometimes one needs to calculate the LD between two SNPs. Select randomly another SNP from your results file and calculate LD with that and the top SNP using LDpair. See how the LD changes when you select a SNP from near/far of your top SNP.

  

You may even create an LD matrix between several SNPs. Select e.g. the ten first SNPs from your results file and calculate an LD matrix between them using LDmatrix.