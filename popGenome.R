## January 2021 Alice Balard Nucleotide diversity for E139, E64, and
## E88 (SNPs per kb) Hyp. of Emanuel: I have the suspicion that the
## original E139 could have been more diverse and therefore responded
## to selection, while the original E64 was less diverse and couldn't
## respond therfore...

## This script calls files that were created based on the code detailed in 
## /SAN/Eimeria_Wild_Genomes/Assembly_Masurca/readme.txt 

library("PopGenome")
setwd("/SAN/Eimeria_Wild_Genomes/Assembly_Masurca")

## E. ferrisi analysis. Let's use E64 as reference genome, and mapps
## E139 pop and E64 pop to it; then get within and between diversity

## Split A VCF File Into Multiple Scaffold-VCFs
## VCF_split_into_scaffolds(VCF.file="vcffiles/E139and64_to_E64gen.vcf", "vcffiles/E64_bothpop64139mapped_scaffolds")

## READ IN DATA
GENOME.class.E139and64_to64gen <- readData("vcffiles/E64_bothpop64139mapped_scaffolds/", format="VCF")

## calculate nucleotide diversity
GENOME.class.E139and64_to64gen <- diversity.stats(GENOME.class.E139and64_to64gen, pi = "True")
mean(GENOME.class.E139and64_to64gen@nuc.diversity.within)
mean(GENOME.class.E139and64_to64gen@nuc.diversity.between)

## Assess populations
get.individuals(GENOME.class.E139and64_to64gen)

GENOME.class.E139and64_to64gen <- set.populations(GENOME.class.E139and64_to64gen,list(c("bamfilesBWA/E139_to_E64gen_coverage.sorted.bam", "bamfilesBWA/E139_to_E64gen_coverage.sorted.bam.2"), c("bamfilesBWA/E64_coverage.sorted.bam", "bamfilesBWA/E64_coverage.sorted.bam.2")))

# Check whether grouping is set correctly
GENOME.class.E139and64_to64gen@region.data@populations2

# Recalculate statistics for populations
GENOME.class.E139and64_to64gen <- diversity.stats(GENOME.class.E139and64_to64gen, pi = "True")

mean(GENOME.class.E139and64_to64gen@nuc.diversity.within)
mean(GENOME.class.E139and64_to64gen@nuc.diversity.between)


GENOME.class <-neutrality.stats(GENOME.class,detail=TRUE)
GENOME.class@Tajima.D
                                        # Each population
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]] 

GENOME.class.E139and64_to64gen@region.stats@nucleotide.diversity

## removing contaminated contigs identified by blobtools

## Load files
E64allmapped <- readLines("../vcffiles/E64allMAPPEDscaffolds.txt")

getScafNum <- function(input){
    data.frame(scaffname = input, scaffnum = 1:length(input))
}

dfE64allmapped = getScafNum(E64allmapped)

## Mapped to apicoplexa
dfE64api = merge(dfE64allmapped,
                  data.frame(scaffname=readLines("../blobtools/E64_blob/scaffolds_Apicomplexa.txt"), all.y=T))

## a lot of well assigned scaffolds can't be found in the vcf file
## e.g. scf7180000052675. Probably no SNPs.
## keep api scaffolds numbers with SNPs
dfE64api = dfE64api[!is.na(dfE64api$scaffnum),]

## Mapped to scaffolds with no hits
dfE64nohit = merge(dfE64allmapped,
                  data.frame(scaffname=readLines("../blobtools/E64_blob/scaffolds_nohit.txt"), all.y=T))

## keep api scaffolds numbers with SNPs
dfE64nohit = dfE64nohit[!is.na(dfE64nohit$scaffnum),]

######### Get scaffolds numbers
scaffoldsAPI64 = dfE64api$scaffnum
scaffoldsAPIorUNKNOWN64 = c(dfE64api$scaffnum, dfE64nohit$scaffnum)    

length(scaffoldsAPI64)
length(scaffoldsAPIorUNKNOWN64)

## calculate nucleotide diversity for each subset
mean(GENOME.class.E139and64_to64gen@nuc.diversity.within)
mean(GENOME.class.E139and64_to64gen@nuc.diversity.within[scaffoldsAPI64])
mean(GENOME.class.E139and64_to64gen@nuc.diversity.within[scaffoldsAPIorUNKNOWN64])

##NB!! The nucleotide diversities have to be devided by GENOME.class@n.sites to give diversities per site.


###### Older analysis (each population studied separately)
## Split A VCF File Into Multiple Scaffold-VCFs
VCF_split_into_scaffolds(VCF.file="../vcffiles/E139.vcf", "../vcffiles/E139_scaffolds")
VCF_split_into_scaffolds(VCF.file="../vcffiles/E88.vcf", "../vcffiles/E88_scaffolds")
VCF_split_into_scaffolds(VCF.file="../vcffiles/E64.vcf", "../vcffiles/E64_scaffolds")

## READ IN DATA
GENOME.class.E139 <- readData("../vcffiles/E139_scaffolds/", format="VCF")
GENOME.class.E88 <- readData("../vcffiles/E88_scaffolds/", format="VCF")
GENOME.class.E64 <- readData("../vcffiles/E64_scaffolds/", format="VCF")

## to get the number of snps per scaffold:
head(GENOME.class.E139@n.biallelic.sites + GENOME.class.E139@n.polyallelic.sites)

GENOME.class.E139@n.biallelic.sites

GENOME.class.E139@n.polyallelic.sites

## calculate nucleotide diversity
GENOME.class.E139 <- diversity.stats(GENOME.class.E139, pi = "True")
mean(GENOME.class.E139@nuc.diversity.within)

GENOME.class.E88 <- diversity.stats(GENOME.class.E88, pi = "True")

mean(GENOME.class.E88@nuc.diversity.within)

GENOME.class.E64 <- diversity.stats(GENOME.class.E64, pi = "True")

mean(GENOME.class.E64@nuc.diversity.within)

############# Part 2: same, but with the contaminated contigs removed from the vcf file

## Load files
E139allmapped <- readLines("../vcffiles/E139allMAPPEDscaffolds.txt")
E88allmapped <- readLines("../vcffiles/E88allMAPPEDscaffolds.txt")
E64allmapped <- readLines("../vcffiles/E64allMAPPEDscaffolds.txt")

getScafNum <- function(input){
    data.frame(scaffname = input, scaffnum = 1:length(input))
}

dfE139allmapped = getScafNum(E139allmapped)
df88allmapped = getScafNum(E88allmapped)
dfE64allmapped = getScafNum(E64allmapped)

## Mapped to apicoplexa
dfE139api = merge(dfE139allmapped,
                  data.frame(scaffname=readLines("../blobtools/E139_blob/scaffolds_Apicomplexa.txt"), all.y=T))
dfE88api = merge(dfE88allmapped,
                  data.frame(scaffname=readLines("../blobtools/E88_blob/scaffolds_Apicomplexa.txt"), all.y=T))
dfE64api = merge(dfE64allmapped,
                  data.frame(scaffname=readLines("../blobtools/E64_blob/scaffolds_Apicomplexa.txt"), all.y=T))

## a lot of well assigned scaffolds can't be found in the vcf file
## e.g. scf7180000052675. Probably no SNPs.
## keep api scaffolds numbers with SNPs
dfE139api = dfE139api[!is.na(dfE139api$scaffnum),]
dfE88api = dfE88api[!is.na(dfE88api$scaffnum),]
dfE64api = dfE64api[!is.na(dfE64api$scaffnum),]

## Mapped to scaffolds with no hits
dfE139nohit = merge(dfE139allmapped,
                    data.frame(scaffname=readLines("../blobtools/E139_blob/scaffolds_nohit.txt"), all.y=T))
dfE88nohit = merge(dfE88allmapped,
                  data.frame(scaffname=readLines("../blobtools/E88_blob/scaffolds_nohit.txt"), all.y=T))
dfE64nohit = merge(dfE64allmapped,
                  data.frame(scaffname=readLines("../blobtools/E64_blob/scaffolds_nohit.txt"), all.y=T))

## keep api scaffolds numbers with SNPs
dfE139nohit = dfE139nohit[!is.na(dfE139nohit$scaffnum),]
dfE88nohit = dfE88nohit[!is.na(dfE88nohit$scaffnum),]
dfE64nohit = dfE64nohit[!is.na(dfE64nohit$scaffnum),]

######### Get scaffolds numbers
scaffoldsAPI139 = dfE139api$scaffnum
scaffoldsAPIorUNKNOWN139 = c(dfE139api$scaffnum, dfE139nohit$scaffnum)    
scaffoldsAPI88 = dfE88api$scaffnum
scaffoldsAPIorUNKNOWN88 = c(dfE88api$scaffnum, dfE88nohit$scaffnum)    
scaffoldsAPI64 = dfE64api$scaffnum
scaffoldsAPIorUNKNOWN64 = c(dfE64api$scaffnum, dfE64nohit$scaffnum)    

length(scaffoldsAPI139)
length(scaffoldsAPIorUNKNOWN139)
length(scaffoldsAPI88)
length(scaffoldsAPIorUNKNOWN88)
length(scaffoldsAPI64)
length(scaffoldsAPIorUNKNOWN64)

## calculate nucleotide diversity for each subset
mean(GENOME.class.E139@nuc.diversity.within)
mean(GENOME.class.E139@nuc.diversity.within[scaffoldsAPI139])
mean(GENOME.class.E139@nuc.diversity.within[scaffoldsAPIorUNKNOWN139])

mean(GENOME.class.E88@nuc.diversity.within)
mean(GENOME.class.E88@nuc.diversity.within[scaffoldsAPI88])
mean(GENOME.class.E88@nuc.diversity.within[scaffoldsAPIorUNKNOWN88])

GENOME.class.E64@nuc.diversity.within

mean(GENOME.class.E64@nuc.diversity.within[scaffoldsAPI64])
mean(GENOME.class.E64@nuc.diversity.within[scaffoldsAPIorUNKNOWN64])
