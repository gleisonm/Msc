
## Tópicos Principais
-  Data Sparse: Remove proteins that are not present in some patients
-  Integração de dados
	- Middle integration
	- [[JDR]]

## Conceitos-Chave
[[Genomic]] / [[Transcriptomic]] / [[Proteomic]] / [[Metabolomic]]    ->          **Fenótipo**         <-    Análise / Experimento / Dados

Ômicas -> 
Integração (Correlação, Similaridade, Network, Bayesian, Fusion, Multivariate ) -> 
Análise (Subtipo de doença, Insights de doença, Predição de biomarcadores) -> 
Resultados

## Resumo
Para otimizar os dados remover dados ausentes dos pacientes ([[Sparse]]) (*Estudar isso*)
Durante a integração precisamos reduzir a dimensionalidade da matriz de cada ômica ([[JDR]])

### Integração
- Early Integration: Integra Dados
- Middle Integration: Integra Dados a partir de um modelo multiômico **Core da Multiômica (Fatores comuns)**
- Late Integration: Integra Resultados

#### Multi-omics Joint Dimentionality Reduction ([[JDR]])
- "PCA da multiômica"
- Cria uma matriz de fatores e amostras a partir de uma matriz específica para cada ômica
- Necessário para clusterização de amostras e análise "bioinformata", com redes, processos, mecanismos...
- Existem diversos tipos de redução de dimensionalidade, ex:
	- JIVE
	- MCIA
	- FA
	- MOFA

##### Joint and Individual Variance Explain ([[JIVE]])
- É um tipo de JDR que diminui o [[chi quadrado]]
- Estrutura de junção + Estruturas Individuais + Ruído
-  Baseado em PCA
-  Quebra a matriz em componentes e encontra fatores em conjuntos
- Encontra entre as matrizes concordâncias e discordâncias e minimizam as diferenças

##### Multiple Co-Inertia Analysis ([[MCIA]])
- Quando você possui mais de um dataset
- Maximiza a covariância

##### Análise de Fatores ([[FA]])
- Maximiza a correlação ou covariância
- Conjuntos ou informações não observáveis (Variáveis Latentes)
- Ferramenta: [[ICluster]]

##### Multi-Omics Factor Analysis ([[MOFA]])
- Parte de um modelo probabilístico Bayesiano

## Links Úteis
- [Benchmarking joint multi-omics dimensionality reduction approaches for the study of cancer](https://www.nature.com/articles/s41467-020-20430-7)
- 

## Tarefas e Ações
- [ ] Revisar Sparse 
- [ ] Revisar JDR

## Dúvidas e Perguntas
-  



#Multiomic #Genomic #Transcriptomic #Proteomic #Metabolomic #class