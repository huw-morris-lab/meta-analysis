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
