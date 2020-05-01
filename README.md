# meta-analysis
Steps to meta-analyse GWAS results from multiple studies.

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

I do this at the same time as getting GWAS results into the right format for FUMA - so there may be some overlap in the steps here in and the GWAS post-processing.

Make combined file if you have run your GWAS on subsets of the data separately, e.g. by chromosome or in subsets of 50k SNPs. 
```
cat GWASresults_* | grep -v 'SNP' > allGWAS_results.txt
```

Read data into R to calculate lambda value, annotate with rsIDs, and export data in correct formats for FUMA and METAL. Note that this script is working with survival (Cox proportional hazard model) GWAS results so the headers of your results file may be different.

Check https://genome.sph.umich.edu/wiki/METAL_Documentation for the correct format of the input files for METAL.

```
R
library(dplyr)
library(data.table)
library(tidyr)

#Read in data
data <- fread("allGWAS_results.txt")

#Label columns - this is from a survival GWAS
colnames(data) <- c("SNP", "Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio", "logrank", "r2")

#Sort by p value
data <- data %>% arrange(Pvalue)

#Check Cox.zphPVal
#This is really, really important to get right. If this p-value is very small, then it means that the variable is time dependent and therefore needs to be treated separately. 
#Check this column if you get any significant hits.
 
#Calculate lambda
p <- data$Pvalue
n <- length(data$Pvalue)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda

#Remove X from the SNP names if you have them - found this was when running GWAS in R.
data <- data %>%
	mutate(SNP = substring(SNP, 2))

#Split SNP name into chr and position and alleles
data_split <- data %>%
	separate(SNP, into = c("chr", "bp", "REF", "ALT", "A1"))
#There are some warnings, I think these are indels
#The ALT allele appears as CN0
#This is fine - in the bim file the A1 is empty so the SNP is a deletion

#Check some instances where the A1 coded by plink --recode is not the same as the ALT column (from VCF)
#We need to use the plink A1 column because this is the allele count that is used for the survival model
#Everything else is the non-effect allele
data_split  <- data_split %>%
  mutate(effect = ifelse(A1 == ALT, ALT,
                          ifelse(A1 == REF, REF, NA)),
         noneffect = ifelse(A1 == ALT, REF,
                          ifelse(A1 == REF, ALT, NA)))


#Remove indels
data_split_noindels <- data_split %>%
	filter(!is.na(effect)) %>%
	filter(!is.na(noneffect))

#Export for FUMA
fwrite(data_split_noindels, "results_FUMA.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Read in frequency file - needed for METAL
freq <- fread("$FILENAME.frq")

freq <- freq %>%
	select(-CHR, -NCHROBS) %>%
	rename(A1_freq = A1,
		A2_freq = A2)

#Merge with results file
data_split_noindels_freq <- data_split_noindels %>%
	mutate(SNP = paste(chr,bp,REF,ALT, sep = ":")) %>%
	left_join(freq, by = "SNP")

#Check that alleles match (there should be no mismatches)
data_split_noindels_freq %>%
	filter(A1!=A1_freq) %>%
	summarise(count = n())

#Read in rsIDs file - we need rsIDs for meta-analysis 
#Because some datasets in our meta-analysis are in hg38 so going to analyse by rsID
rsids <- fread("/data/kronos/mtan/HRC_rs_ids_GRCh37.txt")

#Merge results with rsIDs
data_split_noindels_freq_rsids <- data_split_noindels_freq %>%
	mutate(chrbp = paste(chr,bp, sep =":")) %>%
	inner_join(rsids, by = "chrbp")
#5464833 SNPs that have rsIDs

#Check allele matches
data_split_noindels_freq_rsids <- data_split_noindels_freq_rsids %>%
	mutate(allele_match = ifelse(effect == ALT.y & noneffect == REF.y, "match1",
				ifelse(effect == REF.y & noneffect == ALT.y, "match2", "mismatch")))


#Count allele matches and mismatches
data_split_noindels_freq_rsids %>%
	group_by(allele_match) %>%
	summarise(count = n())

#Filter out allele mismatches
data_split_noindels_freq_rsids_alleleMatch <- data_split_noindels_freq_rsids %>%
	filter(allele_match == "match1" | allele_match == "match2")

#Check that rsIDs are unique
data_split_noindels_freq_rsids_alleleMatch[duplicated(data_split_noindels_freq_rsids_alleleMatch$ID), ]

#There are duplicated elements where the rsIDs are missing
#Remove these
data_split_noindels_freq_rsids_alleleMatch <- data_split_noindels_freq_rsids_alleleMatch %>%
	filter(!is.na(ID)) %>%
	filter(ID!=".")

#Also still duplicated rsIDs where there are positions with multiple alleles e.g. 1:96656801:C:A and 1:96656801:C:T with the same rsID
#Remove these 
#Have ordered data according to p value so should keep one of the pair with the smaller p value
data_split_noindels_freq_rsids_alleleMatch_unique <- data_split_noindels_freq_rsids_alleleMatch %>% distinct(ID, .keep_all = TRUE)

#Format for METAL
export_METAL <- data_split_noindels_freq_rsids_alleleMatch_unique %>%
	select(ID, effect, noneffect, Coeff, se, Pvalue, N, MAF) %>%
	rename(rsid = ID,
		effect_allele = effect,
		noneffect_allele = noneffect,
		beta = Coeff)

#Export for METAL
fwrite(export_METAL, "results_METAL.tab", quote = F, sep = "\t", row.names = F)

q()
n
```

If your file is too large to upload to FUMA, compress it with gzip.
```
gzip results_FUMA.txt
```



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
