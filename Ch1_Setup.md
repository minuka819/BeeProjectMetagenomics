Setting up new Phyloseq Objects
================
Minuka Hewapathirana
2024-07-19

### (1) Loading Packages

A comprehensive list of packages required to perform metagenomics
analysis for *ION torrent* data

``` r
library("knitr")
library("qiime2R") # devtools::install_github("jbisanz/qiime2R")
library("phyloseq")
library("readxl")      
library("tibble")
library("vegan")
library("DESeq2") #BiocManager::install("DESeq2")
library("speedyseq") # remotes::install_github("mikemc/speedyseq") 
library("ape")
library("ggstar")
library("forcats")
library("patchwork")
library("ggpubr")
library("plotROC")
library("viridis")
library("cowplot")
library("ggplot2")
library("microbiome") # BiocManager::install("microbiome")
library("microbiomeutilities")

library("ggtree") # BiocManager::install("ggtree")
library("ggtreeExtra") #install.packages("ggExtra")
library('MicrobiotaProcess') # BiocManager::install("MicrobiotaProcess")
#library("tidytree")
library("file2meco")
library("microeco")
```

### (2) Converting Qiime Artifacts to Phyloseq Objects

Each **Phyloseq** object requires - features (table.qza) - tree
(rooted-tree.qza) - taxonomy (taxonomy.qza) - metadata (“sample.tsv”) :
The metadata file shoudl be in .tsv format

``` r
#Creating a phyloseq object for the ITS1 Reverse Reads 
PS_ITS1R_Global<-qza_to_phyloseq(
  features="D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\ITS2_NEW\\table.qza",
  tree="D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\ITS2_NEW\\rooted-tree.qza",
  taxonomy="D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\ITS2_NEW\\taxonomy.qza",
  metadata = "D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\ITS2_NEW\\ITS2_Global.tsv"
)

#Checking if objects exist
PS_ITS1R_Global
```

phyloseq-class experiment-level object otu_table() OTU Table: \[ 11005
taxa and 251 samples \]: sample_data() Sample Data: \[ 251 samples by 15
sample variables \]: tax_table() Taxonomy Table: \[ 11005 taxa by 7
taxonomic ranks \]: phy_tree() Phylogenetic Tree: \[ 11005 tips and
10923 internal nodes \]: taxa are rows

``` r
# Creating a phyloseq object for the 16S Reverse Reads 
PS_16S_Global<-qza_to_phyloseq(
  features="D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\New_Silva_Global\\table.qza",
  tree="D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\New_Silva_Global\\rooted-tree.qza",
  taxonomy="D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\New_Silva_Global\\taxonomy.qza",
  metadata = "D:\\Grad_School\\Grad_Project\\Qiime_Outputs\\New_Silva_Global\\Global.tsv"
)

PS_16S_Global
```

phyloseq-class experiment-level object otu_table() OTU Table: \[ 11535
taxa and 254 samples \]: sample_data() Sample Data: \[ 254 samples by 14
sample variables \]: tax_table() Taxonomy Table: \[ 11535 taxa by 7
taxonomic ranks \]: phy_tree() Phylogenetic Tree: \[ 11535 tips and
11352 internal nodes \]: taxa are rows

### (3) Obtaining StoneFruit and Blueberry ITS1 Datasets

Goal : Split ITS1 PS Object into two smaller PS objects based on
StoneFruit and Blueberry in BC and ON

- First remove Alberta from ITS1R and 16S objects

- Then split based on Stonefruits and Blueberries respectively

``` r
#removal of AB province
PS_ITS1R_BC_ON <- subset_samples(PS_ITS1R_Global, Province !="AB")

#removal of Wild forage 
PS_ITS1R_BC_ON_SF <- subset_samples(PS_ITS1R_BC_ON, Plant !="Wild forage")
#removal of Blueberry 
PS_ITS1R_BC_ON_SF <- subset_samples(PS_ITS1R_BC_ON_SF, Plant !="Blueberry")
#PS Object with only Stonefruit in BC and ON 
PS_ITS1R_BC_ON_SF
```

phyloseq-class experiment-level object otu_table() OTU Table: \[ 11005
taxa and 138 samples \]: sample_data() Sample Data: \[ 138 samples by 15
sample variables \]: tax_table() Taxonomy Table: \[ 11005 taxa by 7
taxonomic ranks \]: phy_tree() Phylogenetic Tree: \[ 11005 tips and
10923 internal nodes \]: taxa are rows

``` r
#Check metadata file for Stonefruit
metadata_PS_ITS1R_BC_ON_SF <- meta(PS_ITS1R_BC_ON_SF)

kable(head(metadata_PS_ITS1R_BC_ON_SF), digits = 2, align = c(rep("l", 4), rep("c", 4), rep("r", 4)))
```

|     | barcode.seqeunce | linker.primer.seqeunce     | InputFileName              | Description | RunNumber | SampleName  | Number | Province |       Location | Plant | Sample.type | Site | DateCollected | Colony | DNANg |
|:----|:-----------------|:---------------------------|:---------------------------|:------------|:---------:|:-----------:|:------:|:--------:|---------------:|------:|------------:|-----:|:--------------|:-------|:------|
| A01 | CTAAGGTAA        | CGATCTTGGTCATTTAGAGGAAGTAA | A01_9_L001_R1_001.fastq.gz | ITS2_A_B01  |  BCC_R1   | BCCV1-AP-1B |   1    |    BC    | Creston Valley | Apple |       Bread |    1 | Spring 2021   | 1      | 1.16  |
| A02 | TAAGGAGAA        | CGATCTTGGTCATTTAGAGGAAGTAA | A02_9_L001_R1_001.fastq.gz | ITS2_A_B02  |  BCC_R1   | BCCV1-AP-2B |   2    |    BC    | Creston Valley | Apple |       Bread |    1 | Spring 2021   | 2      | 0.74  |
| A03 | AAGAGGATT        | CGATCTTGGTCATTTAGAGGAAGTAA | A03_9_L001_R1_001.fastq.gz | ITS2_A_B03  |  BCC_R1   | BCCV1-AP-3B |   3    |    BC    | Creston Valley | Apple |       Bread |    1 | Spring 2021   | 3      | 0.87  |
| A04 | TACCAAGAT        | CGATCTTGGTCATTTAGAGGAAGTAA | A04_9_L001_R1_001.fastq.gz | ITS2_A_B04  |  BCC_R1   | BCCV2-AP-1B |   4    |    BC    | Creston Valley | Apple |       Bread |    2 | Spring 2021   | 1      | 0.56  |
| A05 | CAGAAGGAA        | CGATCTTGGTCATTTAGAGGAAGTAA | A05_9_L001_R1_001.fastq.gz | ITS2_A_B05  |  BCC_R1   | BCCV2-AP-2B |   5    |    BC    | Creston Valley | Apple |       Bread |    2 | Spring 2021   | 2      | 0.99  |
| A06 | CTGCAAGTT        | CGATCTTGGTCATTTAGAGGAAGTAA | A06_9_L001_R1_001.fastq.gz | ITS2_A_B06  |  BCC_R1   | BCCV2-AP-3B |   6    |    BC    | Creston Valley | Apple |       Bread |    2 | Spring 2021   | 3      | 0.52  |

``` r
#Remove other plant types to obtain Wild Forage and Bluberry subset 
#Cherry, Apple , Peach , Apricot 
PS_ITS1R_BC_ON_BB <-subset_samples(PS_ITS1R_BC_ON, Plant !="Cherry")
PS_ITS1R_BC_ON_BB <-subset_samples(PS_ITS1R_BC_ON_BB, Plant !="Apple")
PS_ITS1R_BC_ON_BB <- subset_samples(PS_ITS1R_BC_ON_BB, Plant !="Peach")
PS_ITS1R_BC_ON_BB <- subset_samples(PS_ITS1R_BC_ON_BB, Plant !="Apricot")

PS_ITS1R_BC_ON_BB
```

phyloseq-class experiment-level object otu_table() OTU Table: \[ 11005
taxa and 82 samples \]: sample_data() Sample Data: \[ 82 samples by 15
sample variables \]: tax_table() Taxonomy Table: \[ 11005 taxa by 7
taxonomic ranks \]: phy_tree() Phylogenetic Tree: \[ 11005 tips and
10923 internal nodes \]: taxa are rows

``` r
#Check metadatafile for blueberry
metadata_ITS1R_BC_ON_BB <- meta(PS_ITS1R_BC_ON_BB)

kable(head(metadata_ITS1R_BC_ON_BB), digits = 2, align = c(rep("l", 4), rep("c", 4), rep("r", 4)))
```

|     | barcode.seqeunce | linker.primer.seqeunce     | InputFileName              | Description | RunNumber | SampleName  | Number | Province |      Location |     Plant | Sample.type | Site | DateCollected | Colony | DNANg |
|:----|:-----------------|:---------------------------|:---------------------------|:------------|:---------:|:-----------:|:------:|:--------:|--------------:|----------:|------------:|-----:|:--------------|:-------|:------|
| A10 | CTGACCGAA        | CGATCTTGGTCATTTAGAGGAAGTAA | A10_9_L001_R1_001.fastq.gz | ITS2_A_B10  |  BCC_R1   | BCFV1-BB-1B |   10   |    BC    | Fraser Valley | Blueberry |       Bread |    1 | 12-May-21     | 1      | 1.83  |
| A11 | TCCTCGAAT        | CGATCTTGGTCATTTAGAGGAAGTAA | A11_9_L001_R1_001.fastq.gz | ITS2_A_B11  |  BCC_R1   | BCFV1-BB-2B |   11   |    BC    | Fraser Valley | Blueberry |       Bread |    1 | 12-May-21     | 2      | 1.36  |
| A12 | TAGGTGGTT        | CGATCTTGGTCATTTAGAGGAAGTAA | A12_9_L001_R1_001.fastq.gz | ITS2_A_B12  |  BCC_R1   | BCFV1-BB-3B |   12   |    BC    | Fraser Valley | Blueberry |       Bread |    1 | 12-May-21     | 3      | 2.51  |
| A13 | TCTAACGGA        | CGATCTTGGTCATTTAGAGGAAGTAA | A13_9_L001_R1_001.fastq.gz | ITS2_A_B13  |  BCC_R1   | BCFV1-BB-4B |   13   |    BC    | Fraser Valley | Blueberry |       Bread |    1 | 12-May-21     | 4      | 4.49  |
| A14 | TTGGAGTGT        | CGATCTTGGTCATTTAGAGGAAGTAA | A14_9_L001_R1_001.fastq.gz | ITS2_A_B14  |  BCC_R1   | BCFV2-BB-1B |   14   |    BC    | Fraser Valley | Blueberry |       Bread |    2 | 12-May-21     | 1      | 5.43  |
| A15 | TCTAGAGGT        | CGATCTTGGTCATTTAGAGGAAGTAA | A15_9_L001_R1_001.fastq.gz | ITS2_A_B15  |  BCC_R1   | BCFV2-BB-2B |   15   |    BC    | Fraser Valley | Blueberry |       Bread |    2 | 12-May-21     | 2      | 3.86  |

### (4) Obtaining StoneFruit and Blueberry 16S Datasets

Goal : Repeat Step 3 for 16S

``` r
#removal of AB province
PS_16S_BC_ON <- subset_samples(PS_16S_Global, Province !="AB")

#meta(PS_16S_BC_ON)

#Stone fruit subset     
PS_16S_BC_ON_SF <- subset_samples(PS_16S_BC_ON, Plant !="Wild forage")
PS_16S_BC_ON_SF <- subset_samples(PS_16S_BC_ON_SF, Plant !="Blueberry")
PS_16S_BC_ON_SF <- subset_samples(PS_16S_BC_ON_SF, Plant !="Kit")

PS_16S_BC_ON_SF
```

phyloseq-class experiment-level object otu_table() OTU Table: \[ 11535
taxa and 140 samples \]: sample_data() Sample Data: \[ 140 samples by 14
sample variables \]: tax_table() Taxonomy Table: \[ 11535 taxa by 7
taxonomic ranks \]: phy_tree() Phylogenetic Tree: \[ 11535 tips and
11352 internal nodes \]: taxa are rows

``` r
#meta(PS_16S_BC_ON_SF)

#Remove these to obtain Wild Forage and Bluberry subset 
# Cherry, Apple , Peach , Apricot 
PS_16S_BC_ON_BB <-subset_samples(PS_16S_BC_ON, Plant !="Cherry")
PS_16S_BC_ON_BB <-subset_samples(PS_16S_BC_ON_BB, Plant !="Apple")
PS_16S_BC_ON_BB <- subset_samples(PS_16S_BC_ON_BB, Plant !="Peach")
PS_16S_BC_ON_BB <- subset_samples(PS_16S_BC_ON_BB, Plant !="Apricot")
PS_16S_BC_ON_BB <- subset_samples(PS_16S_BC_ON_BB, Plant !="Kit")

PS_16S_BC_ON_BB
```

phyloseq-class experiment-level object otu_table() OTU Table: \[ 11535
taxa and 82 samples \]: sample_data() Sample Data: \[ 82 samples by 14
sample variables \]: tax_table() Taxonomy Table: \[ 11535 taxa by 7
taxonomic ranks \]: phy_tree() Phylogenetic Tree: \[ 11535 tips and
11352 internal nodes \]: taxa are rows

``` r
#meta(PS_16S_BC_ON_BB)
```

We have obtained **4 PS objects** that will be used in the following
chapters to perorm metagenomics based visualizations and develop
insights

### *UP NEXT* - Chapter 2 : Metagenomics Analysis of Fungal (ITS1) and bacterial (16S) data to investigate agricultural and invasive pathogen detections in **Stonefruit Systems**
