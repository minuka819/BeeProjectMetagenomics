# Metagenomics Analysis on Bee and Pollen Samples for Stone Fruit and Blueberry Ecosystems in Canada

This project aims at performing analysis on Qiime2 Classified outputs through the Phyloseq Package.
The pipeline is broken down into 3 chapters which can be viewed as github md files.

## Pre-Processing
- The analysis is based on ION torrent data that has been classified through UNITE (fungal) and SILVA (bacterial) via Qiime2
- Output files can be found in the **Output** directory 
- Image files are provided in **Image** directory 

## $${\color{red}Chapter\space 1}$$ : Setup
- Loading packages [Qiime2,Phyloseq,microbiotaprocess,microbiomeutilities...]
- Conversion of Qiime2 .qza files to Phyloseq(PS) object 

## $${\color{red}Chapter\space 2}$$ : Metagenomics Analysis on *Fungal* Communities in Stone Fruit Ecosystems in Canada
- Read distributions
- Rarefaction cureves
- Singleton Processing
- Relative abundance barcharts
- Heatmaps
- Focused Phylogenetic Trees
- Metadata Analysis

## $${\color{red}Chapter\space 3}$$ : Metagenomics Analysis on *Bacterial* Communities in Stone Fruit Ecosystems in Canada 
- Markdown file in development for chapter 3.. coming soon
- Chapter 4 in RMD format is available for bacterial applications on blueberry Ion Torrent data 

## List of Abbreviations
- PS - Phyloseq
- md - markdown
- SF - Stone Fruit
- BB - Blueberry
