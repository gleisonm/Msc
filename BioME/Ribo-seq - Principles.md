## Tópicos Principais
- Várias doenças são afetadas a nível proteico, mas não a nível de transcrito
- Ribossome profiling mede a abundancia dos transcritos ligados aos ribossomos. 
	- É possível identificar mRNA ORFs (de novo) e eventos de transcrição não identificáveis por outros métodos.

## Conceitos-Chave
- [[Transcriptomic]]
- [[Multiômics]]
- [[ORF]] & [[CDS]]
- [[Ribo-seq]]
- Eficiência de sequenciamento e análise:
	- **High Throughput**: NGS - Maior profundidade e cobertura
	- **Low Throughput**: Pequeno conjunto de dados
- A posição relativa aos sítios A e P nos fragmentos protegidos pelo ribossomo pode ser determinada pelas Ribo-seq reads mapeadas pelos start e stop codon.
	- Apenas o sítio P permite a ocupação do **start códon**
	- apenas o sítio A permite a ocupação do **stop códon**
- Tradução:
![[sintese-de-proteina.webp]]

## Resumo
RNA-seq é uma das metodologias mais utilizadas para investigar doenças humanas, entretanto algumas evidencias sugerem que muitas das doenças humanas são **afetadas a nível proteico, mas não a nível de transcrito** , sendo assim algumas análises não podem ser feitas somente por RNA-seq.

A maioria das análises a nível proteico são métodos  **Low Throughput**, já Ribo-seq é um método **High Throughput** que sistematicamente mede a **abundancia dos transcritos ligados aos ribossomos**. 

Atualmente as metodologias de Ribo-seq possuem resolução a nível de nucleotídeos únicos, **sendo possível identificar mRNA ORFs (de novo) e eventos de transcrição não identificáveis por outros métodos.**

### Workflow
A base da Ribo-seq é que **fragmentos de mRNA ligados à ribossomos são protegidos da digestão por nuclease**. Por isso, com utilização de métodos de sequenciamento High Throughput e estratégias computacionais é possível aferir o *ribossome footprint* a partir de transcritos

[[Ribo-seq - keysteps]]
![[Pasted image 20240924115119.jpg]]
>*Schematic illustration of the current workflow of ribosomes. The experiment starts with cell lysis, which isolates and immobilizes the mRNA ribosome complexes, and is followed by nuclease digestion of mRNA sequences that are not protected by associated ribosomes. Purification of the mRNA fragments shielded by the ribosomes is then carried out, followed by standard deep sequencing protocols, such as library preparation.*

RNA-seq mede a quantidade total de mRNA e a regulação transcricional, Ribo-seq, por outro lado, pode ser utilizado para medir os mRNAs que realmente foram traduzidos, sendo mais acurado para medir a expressão proteica do que RNA-seq.

### How to Estimate Protein Abundances by Ribo-seq Data
Tecnicamente, Ribo-seq é utilizado para estimar a abundancia da tradução de mRNAs, mas não o nível de expressão proteico.  [Este estudo estudo](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B16-cells-11-02966) demonstrou que a abundância proteica estimada por Ribo-seq é significativamente correlacionada com a abundância realizada pro espectrofotometria de massa. **However, we should be aware that this correlation only holds true under certain conditions and may not be true for other conditions**
> É Duvidoso porque a espectrofotometria de massa não é capaz de diferenciar proteínas recém sintetizadas com as que já eram existentes, enquanto Ribo-seq mede mRNAs que estão prestes a serem transcritos.

### Computational Methods and Tools
- **Existem dois grandes desafios computacionais:**
	 1. Identificar no fragmento de mRNA protegido pelos ribossomos onde estão os sítios A e P, isso é fundamental para compreender a tradução em escala de códons.
		- Os sítios A P E são separados por 3 nucleotídeos.
		- A posição relativa aos sítios A e P nos fragmentos protegidos pelo ribossomo pode ser determinada pelas Ribo-seq reads mapeadas pelos start e stop codon.
			- Os métodos para identificar esse sítios são baseados no fato de que **apenas o sítio P permite a ocupação do start códon** durante o início da tradução e **apenas o sítio A permite a ocupação do stop códon**.
		- **Problema:** A digestão pela RNase é estocástica, então o tamanho do fragmento protegido pelo ribossomo pode variar, por isso é tecnicamente difícil de determinar os sítios A e P de amostars com diferentes tamanhos de fragmento.
			- Para contornar esse problema um estudo utilizou o fato de que o sítio A nos fragmentos protegidos por ribossomo devem residir dentro da região CDS para criar uma função de maximiza o número de fragmentos com o sítio A caindo na região CDS.
	2. **Ribo-seq data normalization density**: A read de Ribo-seq em um transcrito de referência (.gtf) é caracterizado por gaps no alinhamento, picos de alta densidade devido a artefatos técnicos e pausas de ribossomos.
		- O método computacional simples e robusto é o RUST, que reduz o bias de densidade para que a normalização seja feita de forma correta

##### Quantification of Ribossome-bounded transcripts:

| Methods and Website                                                                                                                                                                                                                                                                                                                                                                          | Features                                                                                                                                    |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| riboSeqR: An R/Bioconductor package that provides a set of programs for processing and visualization of Ribo-seq data.  <br>[https://ribogalaxy.genomicsdatascience.ie/](https://ribogalaxy.genomicsdatascience.ie/) [[20](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B20-cells-11-02966)]                                                                                        | Provides visualization of data at sub-codon resolution in the context of single transcripts.                                                |
| Plastid: A user-friendly, generalized analytical pipeline tool that enables users to manipulate data nucleotide by nucleotide robustly and easily and that is not limited to specific experimental regimes or analytical workflows.  <br>[https://plastid.readthedocs.io](https://plastid.readthedocs.io/) [[21](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B21-cells-11-02966)]. | Extensibility and flexibility across assays while remaining user friendly.                                                                  |
| RUST: A smoothing transformation-based approach for Ribo-seq normalization in the presence of heterogeneous noise.  <br>[https://lapti.ucc.ie/rust/](https://lapti.ucc.ie/rust/) [[22](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B22-cells-11-02966)].                                                                                                                           | Performs better in presence of sporadic heterogeneous noise than the previous methods.                                                      |
| mQC: A tool for visualizing quality and data exploration after mapping.  <br>[https://github.com/Biobix/mQC](https://github.com/Biobix/mQC) [[23](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B23-cells-11-02966)].                                                                                                                                                                | Applies the P site offsets before plotting to inspect ribosomal framing and triplet periodicity more elaborately than other existing tools. |
| GWIPS-viz: An online genome browser for checking quality measures or discovering authentic new information from ribosome profiling data.  <br>[https://gwips.ucc.ie/](https://gwips.ucc.ie/) [[24](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B24-cells-11-02966)].                                                                                                               | A Ribo-seq genome browser for data visualization                                                                                            |
| RiboVIEW: A computational pipeline for visualization, quality control, and statistical analysis of ribosome profiling data.  <br>[https://github.com/carinelegrand/RiboVIEW](https://github.com/carinelegrand/RiboVIEW) [[25](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B25-cells-11-02966)].                                                                                    | Focuses on checking quality measures.                                                                                                       |
| Trips-Viz: A graphical tools for exploring properties of collection of ORFs.  <br>[https://trips.ucc.ie/](https://trips.ucc.ie/) [[26](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B26-cells-11-02966)].                                                                                                                                                                           | Provides visualization of data at sub codon resolution in the context of single transcripts.                                                |
##### Differential translation analysis & Identification of A and P site location

| Methods and Website                                                                                                                                                                                                                                                                                                                                                             | Features                                                                                                                  |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Riborex: A linear model-based tool for identification of differential translation from Ribo-seq data.  <br>[https://github.com/smithlabcode/riborex](https://github.com/smithlabcode/riborex) [[35](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B35-cells-11-02966)].                                                                                                 | Faster than all existing methods and employs robust software implementations for the underlying statistical calculations. |
| Anota: An R/Bioconductor package that implements analysis of partial variance (APV) to identify differential translation.  <br>[https://github.com/ChrOertlin/anota2seq/](https://github.com/ChrOertlin/anota2seq/) [[36](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B36-cells-11-02966)].                                                                           | Using APV instead of log ratio approach for detecting translation changes.                                                |
| Babel: An errors-in-variables regression model-based framework to compare ribosome associations within and between conditions based on an errors-in-variables regression model.  <br>[https://github.com/olshena/babel](https://github.com/olshena/babel) [[37](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B37-cells-11-02966)].                                     | Model is more flexible and combines P-values across independent tests.                                                    |
| RiboDiff: A linear model-based framework for detecting changes of mRNA translation efficiency across experimental conditions.  <br>[http://bioweb.me/ribodiff](http://bioweb.me/ribodiff)  <br>[http://github.com/ratschlab/ribodiff](http://github.com/ratschlab/ribodiff) [[38](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B38-cells-11-02966)].                   | Facilitating comparisons of RF abundance by taking mRNA abundance variability as a confounding factor.                    |
| Xtail: An analysis pipeline to detect differentially translated genes.  <br>[https://github.com/xryanglab/xtail](https://github.com/xryanglab/xtail) [[39](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B39-cells-11-02966)].                                                                                                                                          | A more sophisticated method for domination on limitations, such as high-false discoveries and low sensitivities.          |
| RiboProfiling: An R/Bioconductor package that provides a full pipeline to cover all key steps for facilitating the analysis of Ribo-seq experiments and ribosome footprints.  <br>[https://github.com/alenzhao/RiboProfiling](https://github.com/alenzhao/RiboProfiling)  <br>[[40](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B40-cells-11-02966)].                 | Utilizes multiple R packages to handle datasets easily.                                                                   |
| RiboA: A user-friendly web application that identifies A site locations and generates read density profiles.  <br>Website: [https://a-site.vmhost.psu.edu/](https://a-site.vmhost.psu.edu/)  <br>[https://github.com/obrien-lab/aip_web_docker](https://github.com/obrien-lab/aip_web_docker) [[41](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B41-cells-11-02966)]. | The most accurate identifier compared to other tools.                                                                     |
| riboWaltz: An R package for the identification of the ribosome P site, analysis, and visual inspection of ribosome profiling data.  <br>[https://github.com/LabTranslationalArchitectomics/RiboWaltz](https://github.com/LabTranslationalArchitectomics/RiboWaltz) [[42](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B42-cells-11-02966)].                            | Addresses issue of time limitation and data preprocessing.                                                                |
| RiboToolkit: A freely available, web-based service to centralize Ribo-seq data analyses, codon occupancy, and translation efficiency analysis.  <br>[http://rnainformatics.org.cn/RiboToolkit/](http://rnainformatics.org.cn/RiboToolkit/) [[43](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B43-cells-11-02966)].                                                    | Addresses the lacking integrated tool and easy-to-use integrated tool to analyze Ribo-seq data.                           |
| RiboTools: An open-source Galaxy tool used to evaluate codon occupancy at a specific ribosome site and for translation readthrough events.  <br>[https://testtoolshed.g2.bx.psu.edu/view/rlegendre/ribo_tools](https://testtoolshed.g2.bx.psu.edu/view/rlegendre/ribo_tools) [[44](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/#B44-cells-11-02966)].                  | Facilitates complete qualitative analysis.                                                                                |

## Dúvidas e Perguntas
- Qual o tamanho do fragmento que fica dentro do Ribossomo?
	- O ribossomo possui 3 sítios que guardam o mRNA, todos esses sítios estão sendo sequenciados?

## Links Úteis
- [Tracing Translational Footprint by Ribo-Seq: Principle, Workflow, and Applications to Understand the Mechanism of Human Diseases](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/)
- [Link relacionado 2](#)

## Tarefas e Ações
- [ ]  

#riboseq #article #Transcriptomic