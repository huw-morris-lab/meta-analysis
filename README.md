# meta-analysis
Meta-analysis steps for GWAS

Date: May 2020

Last updated: 01/05/2020

Authors: Manuela Tan

## General description and purpose

Meta-analysis for GWAS results

This covers:
1. Getting your data in the right format 
2. Meta-analysis in METAL
3. Filtering meta-analysis results and getting data in format for FUMA

Read the METAL guide https://genome.sph.umich.edu/wiki/METAL_Documentation


# 1. Format data for METAL


# 2. Run meta-analysis in METAL

If you have not already, download the METAL software from http://csg.sph.umich.edu/abecasis/Metal/download/

Downloaded to your local kronos folder. Uuncompress the gz file.
```
tar xvzf Linux-metal.tar.gz 
```

Then make your METAL script. Run your script using
```
/data/kronos/mtan/software/metal/generic-metal/metal metaanalysis_script.txt
```

SCRIPT - following format recommended https://genome.sph.umich.edu/wiki/METAL_Documentation.
Saved as metaanalysis_script.txt
```
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N 

# LOAD INPUT FILES

# Enable Genomic control correction
GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER rsid
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS STUDY1.tab

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER rsid
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS STUDY2.tab

# === DESCRIBE AND PROCESS THE THIRD INPUT FILE ===
MARKER rsid
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS STUDY3.tab

# === DESCRIBE AND PROCESS THE FOURTH INPUT FILE ===
MARKER rsid
ALLELE effect_allele noneffect_allele
FREQ   MAF
EFFECT beta
STDERR se
PVALUE Pvalue
WEIGHT N 
PROCESS STUDY4.tab

# LOAD HOWEVER MANY INPUT FILES YOU HAVE

OUTFILE OUTPUT_META .tbl
ANALYZE HETEROGENEITY

QUIT
```

# 3. Post-analysis filtering

I have followed standard filtering criteria
e.g. Iwaki PD Progression GWAS https://onlinelibrary.wiley.com/doi/full/10.1002/mds.27845


```
R
library(data.table)
library(dplyr)

#Read in meta analysis results
data <- fread("OUTPUT_META.tbl")

#Filter out SNPs which are not in all the studies
#This is up to you whether you set a minimum N of patients that have been analysed
#I filtered by the number of studies 
data_filtered <- data %>%
	filter(HetDf == 3)

#Sort by p value
data_filtered_sorted <- data_filtered %>%
	arrange(`P-value`)

#Filter out SNPs with HetPVal < 0.05 (Cochran's Q-test for heterogeneity)
#Also filter out SNPs with HetISq > 80
data_filtered_sorted_het <-data_filtered_sorted %>%
	filter(HetPVal > 0.05) %>%
	filter(HetISq < 80)

#Check MAF variability - remove variants with MAF variability > 15%
data_filtered_sorted_het_MAF <- data_filtered_sorted_het %>%
	mutate(MAF_variability = MaxFreq - MinFreq) %>%
	filter(MAF_variability <= 0.15)

#Export for FUMA
export_FUMA <- data_filtered_sorted_het_MAF %>%
	select(MarkerName, `P-value`, Allele1, Allele2, Effect, StdErr, TotalSampleSize) %>%
	rename(rsID = MarkerName,
		pval = `P-value`)
		
fwrite(export_FUMA, "metaanalysis_FUMA.txt", quote = F, row.names = F, col.names = T, sep = "\t")
```


Then upload to FUMA to plot results.
