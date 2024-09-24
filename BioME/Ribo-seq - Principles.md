## Tópicos Principais
- Várias doenças são afetadas a nível proteico, mas não a nível de transcrito
- Ribossome profiling mede a abundancia dos transcritos ligados aos ribossomos. 
	- É possível identificar mRNA ORFs (de novo) e eventos de transcrição não identificáveis por outros métodos.

## Conceitos-Chave
- [[Transcriptomic]]
- [[Ribo-seq]]
- Eficiência de coleta e análise:
	- **High Throughput**: NGS - Maior profundidade e cobertura
	- **Low Throughput**: Pequeno conjunto de dados
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



## Dúvidas e Perguntas
- Qual o tamanho do fragmento que fica dentro do Ribossomo?
	- O ribossomo possui 3 "poços" que guardam o mRNA, todos esses poços estão sendo sequenciados?

## Links Úteis
- [Tracing Translational Footprint by Ribo-Seq: Principle, Workflow, and Applications to Understand the Mechanism of Human Diseases](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9562884/)
- [Link relacionado 2](#)

## Tarefas e Ações
- [ ] Revisar 

#riboseq #article #Transcriptomic