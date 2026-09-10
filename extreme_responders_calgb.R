## MSK exceptional responders CG vs GC_Bev in metastatic urothelial cancer -- CALGB 90601

library(tidyverse)
library(data.table)
library(here)
library(wesanderson)
library(ggpubr)

calgb_maf <- fread((here('data','mut_somatic.maf.txt')))
calbg_clinic <- fread((here('data','clinical_info.csv')))

## exploration of the data
str(calgb_maf) # to identify responders vs non-responders, Sample_ID in calgb_clinic == Tumor_Sample_Barcode in calgb_maf
unique(calgb_maf$Tumor_Sample_Barcode) # 22 samples

# clinical files
str(calbg_clinic)
calbg_clinic$tmc <- as.numeric(calbg_clinic$tmc)
calbg_clinic$mutation_count <- as.numeric(calbg_clinic$mutation_count)
calbg_clinic$sample_coverage <- as.numeric(calbg_clinic$sample_coverage)

responders <- calbg_clinic[calbg_clinic$outcome == 'response', ] 
resp_id <- pull(responders, sample_id)
resp_maf <- calgb_maf %>% filter(Tumor_Sample_Barcode %in% resp_id)

resistents <- calbg_clinic[calbg_clinic$outcome == 'resistance', ] 
resist_id <- pull(resistents, sample_id)
resist_maf <- calgb_maf %>% filter(Tumor_Sample_Barcode %in% resist_id)

# loading other genomic results 
alignment_qc <- fread((here('data','alignment_qc.txt')))
cna_armlevel <- fread((here('data','cna_armlevel.txt')))
cna_facets_run_info.txt <- fread((here('data','cna_facets_run_info.txt'))) ## has purity and ploidy
cna_genelevel.txt <- fread((here('data','cna_genelevel.txt'))) # segs
cna_hisens_run_segmentation.seg.txt <- fread((here('data','cna_hisens_run_segmentation.seg.txt'))) # loc start, loc end, seg mean
cna_purity_run_segmentation.seg.txt <- fread((here('data','cna_purity_run_segmentation.seg.txt'))) # same as above
concordance_qc.txt <- fread((here('data','concordance_qc.txt')))
contamination_qc.txt <- fread((here('data','contamination_qc.txt'))) # highest contamination 1.4
DNA.IntegerCPN_CI.txt <- fread((here('data','DNA.IntegerCPN_CI.txt'))) # mismatchseq1 and 2
HLAlossPrediction_CI.txt <- fread((here('data','HLAlossPrediction_CI.txt'))) 

sample_data.txt <- fread((here('data','sample_data.txt'))) # for all samples purity, ploidy, WGD, MSI scores, SBS with p values
write.table(sample_data.txt , file= "sample_data.csv", sep = ",", row.names = FALSE, quote = FALSE)

## Calculate and plot TMB for Extreme responders vs extreme resistance
variants_tmb <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation","Translation_Start_Site",
                  "Nonstop_Mutation")
variants_silent <- c("Silent") 

tmb_resp <- resp_maf%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)

tmb_resist<-  resist_maf%>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        mutate(tmb= n/30)%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n, tmb)

## for boxplots
tmb_resp$Tumor_Sample_Barcode[c(1:10)] <- "Exceptional Responders" ## sounld be 11 samples?? --> s_C_J3LWLT_P001_d is missing from the mafs df. Also not present in the portal CCS_V2Y3KK1M/somatic/ or /bam folder
tmb_resist$Tumor_Sample_Barcode[c(1:12)] <- "Non-Responders"

## TMB for tumor only cohorts
tmb_two_cohorts <- bind_rows(tmb_resp, tmb_resist)
tmb_two_cohorts$Tumor_Sample_Barcode <- factor(tmb_two_cohorts$Tumor_Sample_Barcode, levels = c("Exceptional Responders", "Non-Responders"))

# plot TMB responders vs non-responders
cbPalette <- c("#D55E00", "#0072B2")
plot_tmb <- tmb_two_cohorts %>% ggplot(aes(x= Tumor_Sample_Barcode, y = tmb)) +
        geom_boxplot(aes(fill=Tumor_Sample_Barcode)) +
        scale_fill_manual(values = cbPalette) +
        xlab("Cohorts") +
        ylab("TMB / Mb") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20),  axis.text.x=element_blank(), axis.ticks.x=element_blank())  
plot_tmb

wilcox.test(tmb~Tumor_Sample_Barcode, data = tmb_two_cohorts) # W = 79.5, p-value = 0.21


## Calculate TMB for comut

calgb_silent_synonymous <- calgb_maf %>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_silent)%>%
        filter(`Variant_Classification %in% variants_silent`== "TRUE") %>%
        dplyr::select(Tumor_Sample_Barcode, silent_muts = n) %>%
        dplyr::rename(Synonymous = silent_muts)

calgb_nonsynonymous <- calgb_maf %>%
        group_by(Tumor_Sample_Barcode)%>%
        dplyr::count(Variant_Classification %in% variants_tmb)%>%
        filter(`Variant_Classification %in% variants_tmb`== "TRUE")%>%
        dplyr::select(Tumor_Sample_Barcode, non_silent_muts = n) %>%
        dplyr::rename(Nonsynonymous = non_silent_muts)

calgb_silent_nonsil <- merge(calgb_nonsynonymous, calgb_silent_synonymous, by= "Tumor_Sample_Barcode")
#write.table(calgb_silent_nonsil, file= "calgb_syn_nonsyn.csv", sep = ",", row.names = FALSE, quote = FALSE)
################################################################################ ==================
## Plot fraction of genome altered in exceptional responders vs extreme resistance
fga_resp <- responders %>% select(sample_id, fraction_genome_altered)
fga_resist <- resistents %>% select(sample_id, fraction_genome_altered)

## for boxplots
fga_resp$sample_id[c(1:10)] <- "Exceptional Responders" 
fga_resist$sample_id[c(1:12)] <- "Non-Responders"

## Fraction of genome altered in the two cohorts
fga_two_cohorts <- bind_rows(fga_resp, fga_resist)
fga_two_cohorts$sample_id <- factor(fga_two_cohorts$sample_id, levels = c("Exceptional Responders", "Non-Responders"))
fga_two_cohorts
# plot fraction of genome altered
plot_fga <- fga_two_cohorts %>% ggplot(aes(x= sample_id, y = fraction_genome_altered)) +
        geom_boxplot(aes(fill=sample_id)) +
        scale_fill_manual(values = cbPalette) +
        xlab("Cohorts") +
        ylab("Percent Genome Altered") +
        theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 16), panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 20), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 20),  axis.text.x=element_blank(), axis.ticks.x=element_blank())  
plot_fga

wilcox.test(fraction_genome_altered~sample_id, data = fga_two_cohorts) # W = 58, p-value = 0.6505


ggarrange(plot_tmb, NULL,  plot_fga,
          nrow = 1, widths = c(1, 0.5, 1), 
          legend = "right", common.legend = TRUE )
################################################################################ ==================

nrow(calgb_maf) # 7348 mutations

#A Remove noncoding variants and keep Missense, Nonsense, Splice Site, FrameShift indels, indels and nonsop mutations
final_maf <- calgb_maf %>% filter(Variant_Classification != "Intron", Variant_Classification != "Silent", Variant_Classification != "intron", Variant_Classification != "3'UTR", Variant_Classification != "5'UTR", Variant_Classification != "IGR", Variant_Classification != "RNA", Variant_Classification != "5'Flank", Variant_Classification != "5'utr", Variant_Classification != "Start_Codon_Del", Variant_Classification != "igr", Variant_Classification != "lincRNA", Variant_Classification != "3'utr", Variant_Classification != "Start_Codon_SNP", Variant_Classification != "De_novo_Start_InFrame", Variant_Classification != "De_novo_Start_OutOfFrame", Variant_Classification != "5'flank") #cuts from 7348 to 5274 mutations

write.table(final_maf, file= "mergedMaf.tsv", sep = "/t", row.names = FALSE, quote = FALSE)

calbg_clinic <- rename(calbg_clinic, Tumor_Sample_Barcode = sample_id)
final_maf_clinic <- left_join(final_maf, calbg_clinic, by = 'Tumor_Sample_Barcode')
str(final_maf_clinic$outcome)


#Count number of samples with putative pathogenic mutations in genes including types "Missense, Nonsense, Splice Site, FrameShift indels, indels and nonsop mutations"
count_mutations <- final_maf_clinic %>% dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, outcome) %>% unique() %>% group_by(Hugo_Symbol, outcome) %>% summarise(count = n()) %>% arrange(desc(count))
count_mutations_spread <- spread(count_mutations, outcome, count) %>% mutate(resistance = replace_na(resistance, 0)) %>% mutate(response = replace_na(response, 0)) %>% arrange(desc(response))


# Fisher's exact test on number of alterations per gene

# Initialize function to create 2x2 contingency matrix
make_and_run_fisher <- function(R_mut, NR_mut, total_R, total_NR){
        # create matrix
        temp_fisher_matrix <- matrix(c(R_mut, total_R - R_mut, NR_mut, total_NR - NR_mut), nrow = 2, ncol = 2)
        
        # run fisher
        temp_p_value <- fisher.test(temp_fisher_matrix)$p.value
        return(temp_p_value)
}
## Run fisher's exact on all genes
#spread_response_gene_count
post_fisher_df <- count_mutations_spread %>% rowwise() %>% 
        mutate(fisher_p_result = make_and_run_fisher(response, resistance, total_R = 11, total_NR = 12)) %>% arrange(fisher_p_result)
post_fisher_df <- post_fisher_df %>% 
        mutate(n_patients = resistance + response) %>% 
        mutate(Enriched_in = ifelse(response > resistance, 'Exceptional Responders', 'Non-Responders')) %>% 
        arrange(fisher_p_result)

## Correct for multiple hypothesis testing
FDR_corrected_df <- cbind(post_fisher_df, q_val = p.adjust(post_fisher_df$fisher_p_result, "fdr")) %>% 
        mutate(n_patients = resistance + response) %>% 
        mutate(Enriched_in = ifelse(response > resistance, 'Exceptional Responders', 'Non-Responders')) %>% 
        arrange(q_val)

#No genes significant by multiple hypothesis testing across 4039 genes...really only 2 genes significant by Fishers at P-value < 0.05, FRY and MGA

# https://www.publichealth.columbia.edu/research/population-health-methods/false-discovery-rate

## Prepare files for Mutsig --- https://www.biostars.org/p/164608/
maf_mutsig <- final_maf %>% select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, Variant_Classification, Tumor_Sample_Barcode)

write.table(maf_mutsig , file= "maf_mutsig.tsv", sep="\t", row.names = FALSE, quote = FALSE)

## Ran MutSig in GenePattern, v 1.3 https://cloud.genepattern.org/gp/pages/index.jsf?jobid=395853&openVisualizers=true&openNewWindow=false

#Plot Fisher's exact Pvalues versus MutSigCV2 output 
calgb_maf <- fread((here('data','mut_somatic.maf.txt')))

MutSig_genes <- fread((here('data','mut_somatic.sig_genes.txt')))
colnames(MutSig_genes)[1] <- "Hugo_Symbol"
colnames(MutSig_genes)[14] <- "mutsig_p"

MutSig_pvals <- merge(select(MutSig_genes, Hugo_Symbol, mutsig_p), select(post_fisher_df, Hugo_Symbol, fisher_p_result, n_patients, Enriched_in))

MutSig_pvals

#Used p-threshold of 0.05 in order to color FRY and MGA
ggplot() +
        geom_point(data = MutSig_pvals, aes(x = -log10(mutsig_p), y = -log10(fisher_p_result), size = n_patients)) +
        geom_point(data = MutSig_pvals[MutSig_pvals$fisher_p_result < 0.05, ],
                   aes(x = -log10(mutsig_p), y = -log10(fisher_p_result), color = Enriched_in, size = n_patients)) +
       scale_color_manual(values = c("darkred", "blue"), name = "Enriched in:") + # cbPalette
        scale_size_continuous(name = 'Patients with alteration') +
        geom_text(data = MutSig_pvals[MutSig_pvals$fisher_p_result < 0.05, ],
                  aes(x = -log10(mutsig_p), y = -log10(fisher_p_result), label = Hugo_Symbol), position = position_nudge(x = -1.8)) +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') + theme_bw() +
        xlab('MutSigCV significance [-log10(p-val)]') + ylab('Responder significance [-log10(p-val)]') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

######## ========================================================================================================================================== #########
## Compare fraction of patients with FRY and MGA mutations across TGCA and Faltas Nat Gen 2016 cohorts

# TCGA
setwd("/Users/filipecarvalho/VAL/TCGA/") ## https://www.biostars.org/p/313063/ to interpret TCGA barcodes
tcga_maf <- read.csv("mc3.blca.maf", sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))

final_tcga_maf <- tcga_maf %>% filter(Variant_Classification != "Intron", Variant_Classification != "Silent", Variant_Classification != "intron", Variant_Classification != "3'UTR", Variant_Classification != "5'UTR", Variant_Classification != "IGR", Variant_Classification != "RNA", Variant_Classification != "5'Flank", Variant_Classification != "5'utr", Variant_Classification != "Start_Codon_Del", Variant_Classification != "igr", Variant_Classification != "lincRNA", Variant_Classification != "3'utr", Variant_Classification != "Start_Codon_SNP", Variant_Classification != "De_novo_Start_InFrame", Variant_Classification != "De_novo_Start_OutOfFrame", Variant_Classification != "5'flank") # cuts from 155232 to 102601  

count(final_maf[final_maf$Hugo_Symbol== c("FRY"), ]) # 6 mutations, 5 patients
count(final_tcga_maf[final_tcga_maf$Hugo_Symbol== c("FRY"), ]) #33


### Faltas Nat Gen 2016 cohort
setwd("~/Dropbox (Partners HealthCare)/msk_calgb/data/blca_cornell_2016/")
faltas_maf <- read.csv("data_mutations.txt", sep = "\t",header = TRUE, comment.char = "#", skipNul = FALSE, fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE, na.strings = c("","NA"))

final_faltas_maf <- faltas_maf %>% filter(Variant_Classification != "Intron", Variant_Classification != "Silent", Variant_Classification != "intron", Variant_Classification != "3'UTR", Variant_Classification != "5'UTR", Variant_Classification != "IGR", Variant_Classification != "RNA", Variant_Classification != "5'Flank", Variant_Classification != "5'utr", Variant_Classification != "Start_Codon_Del", Variant_Classification != "igr", Variant_Classification != "lincRNA", Variant_Classification != "3'utr", Variant_Classification != "Start_Codon_SNP", Variant_Classification != "De_novo_Start_InFrame", Variant_Classification != "De_novo_Start_OutOfFrame", Variant_Classification != "5'flank") ## cuts SNVs from 19392 to 13495

count(final_maf[final_maf$Hugo_Symbol== c("FRY"), ]) # 6
count(final_faltas_maf[final_faltas_maf$Hugo_Symbol== c("FRY"), ]) # 0
count(final_faltas_maf[final_faltas_maf$Hugo_Symbol== c("MGA"), ]) # 2
count(final_faltas_maf[final_faltas_maf$Hugo_Symbol== c("TP53"), ]) 
count(final_maf[final_maf$Hugo_Symbol== c("TP53"), ]) 

unique(final_faltas_maf$Tumor_Sample_Barcode) #72

resp_genes <- c('FRY')

test_df <- final_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% resp_genes)%>%
        filter(`Hugo_Symbol %in% resp_genes`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 5 patients, all in responders
test_df$Tumor_Sample_Barcode[c(1:5)] <- "CALGB" 
sum(test_df$patients)/10

test_faltas <- final_faltas_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% resp_genes)%>%
        filter(`Hugo_Symbol %in% resp_genes`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # Zero patients!

test_tcga <- final_tcga_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% resp_genes)%>%
        filter(`Hugo_Symbol %in% resp_genes`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 32 patients
sum(test_tcga$patients)/411 # 0.07785888

test_FRY <- data.frame(cohort= c("CALGB", "TCGA", "Faltas et al."),
                     percent_patients=c(50, 7.8, 0))

ggplot(test_FRY, aes(x= cohort, y= percent_patients, fill= cohort )) +
        geom_bar(stat="identity") + # https://www.learnbyexample.org/r-bar-plot-ggplot2/
        ggtitle("FRY") +
        ylab("Patients with somatic FRY mutations (%)") + xlab("Cohort") +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,50,0), size = 24), legend.title = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 18), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 18)) +
        scale_fill_manual(values = c('#990000', '#009966', '#009966')) # http://www.visibone.com/ for color coding
        
test_fisher_FRY <- matrix(c(5,5,72,0),nrow =2)
fisher.test(test_fisher_FRY) #p-value = 9.236e-06

test_fisher_FRY_tcga <- matrix(c(5,5,379,32),nrow =2)
fisher.test(test_fisher_FRY_tcga) # p-value = 0.0007329

## test TP53 as control
p53 <-  c('TP53')

p53_calgb <- final_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% p53)%>%
        filter(`Hugo_Symbol %in% p53`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 11 patients total, 5 in responders --> 50% responders

p53_faltas <- final_faltas_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% p53)%>%
        filter(`Hugo_Symbol %in% p53`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 41 patients 
sum(p53_faltas$patients)/72 # 0.5694444

p53_tcga <- final_tcga_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% p53)%>%
        filter(`Hugo_Symbol %in% p53`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 201
sum(p53_tcga$patients)/411 # 0.4890511

test_TP53 <- data.frame(cohort= c("CALGB", "TCGA", "Faltas et al."),
                       percent_patients=c(50, 49, 57))

ggplot(test_TP53, aes(x= cohort, y= percent_patients, fill= cohort )) +
        geom_bar(stat="identity") + 
        ggtitle("TP53") +
        ylab("Patients with somatic TP53 mutations (%)") + xlab("Cohort") +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,50,0), size = 24), legend.title = element_blank(), legend.position = 'none', panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 18), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 18)) +
        scale_fill_manual(values = c('#990000', '#009966', '#009966'))


test_fisher_TP53 <- matrix(c(5,5,31,41),nrow =2)
fisher.test(test_fisher_TP53) #p-value = 0.7421

test_fisher_TP53_tcga <- matrix(c(5,5,210,201),nrow =2)
fisher.test(test_fisher_TP53_tcga) # p-value = 1

## test MGA as a potential target
mga <-  c('MGA')

mga_calgb <- final_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% mga)%>%
        filter(`Hugo_Symbol %in% mga`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 4 patients total, all responders --> 40% responders

mga_faltas <- final_faltas_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% mga)%>%
        filter(`Hugo_Symbol %in% mga`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 2 patients 
sum(mga_faltas$patients)/72 # 0.02777778

mga_tcga <- final_tcga_maf %>%
        group_by(Tumor_Sample_Barcode, Hugo_Symbol)%>%
        count(Hugo_Symbol%in% mga)%>%
        filter(`Hugo_Symbol %in% mga`== "TRUE")%>%
        count(Tumor_Sample_Barcode) %>% 
        select(Tumor_Sample_Barcode, Hugo_Symbol,  patients = n) %>% 
        distinct() # 25
sum(mga_tcga$patients)/411 # 0.06

test_MGA <- data.frame(cohort= c("CALGB", "TCGA", "Faltas et al."),
                        percent_patients=c(40, 6, 2.8))

ggplot(test_MGA, aes(x= cohort, y= percent_patients, fill= cohort )) +
        geom_bar(stat="identity") + 
        ggtitle("MGA") +
        ylab("Patients with somatic MGA mutations (%)") + xlab("Cohort") +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,50,0), size = 24), legend.title = element_blank(), legend.position = 'none', panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black", size = 18), axis.title.x = element_blank(), axis.title.y = element_text(colour="black", size = 18)) +
        scale_fill_manual(values = c('#990000', '#009966', '#009966'))


test_fisher_MGA <- matrix(c(6,4,70,2),nrow =2)
fisher.test(test_fisher_MGA) # p-value = 0.001585

test_fisher_MGA_tcga <- matrix(c(6,4,386,25),nrow =2)
fisher.test(test_fisher_MGA_tcga) # p-value = 0.002879


## Look into BAM in IGV to make sure the FRY are not artifacts 
fry_mut <- calgb_maf %>% filter(Hugo_Symbol == 'FRY') ## 7 mutations
fry_mut_final <- final_maf %>% filter(Hugo_Symbol == 'FRY') ## Cuts off one mutation, the silent mutation in K8F2FC

## Workflwo to load BAMs in IGV is on a separate .txt file "igv_notes_calgb_jan22"

###### =========================================================================================================================================== ######

## CNA with FACETS

# responder
#armlevel_s_C_4ETY8N_P001_d <- fread((here('data','s_C_4ETY8N_P001_d__s_C_4ETY8N_N001_d.arm_level.txt')))
genelevel_s_C_4ETY8N_P001_d <- fread((here('data','s_C_4ETY8N_P001_d__s_C_4ETY8N_N001_d.gene_level.txt'))) ## the way to go to find CNA in specific genes!
genelevel_s_C_4ETY8N_P001_d[genelevel_s_C_4ETY8N_P001_d$cn_state != 'DIPLOID',]
genelevel_s_C_4ETY8N_P001_d[genelevel_s_C_4ETY8N_P001_d$gene == 'CDKN2A',]
genelevel_s_C_4ETY8N_P001_d[genelevel_s_C_4ETY8N_P001_d$cn_state == 'DIPLOID',]

# non responder
genelevel_s_C_90MEU5_P001_d <- fread((here('data', 's_C_90MEU5_P001_d__s_C_90MEU5_N001_d.gene_level.txt')))
genelevel_s_C_90MEU5_P001_d[genelevel_s_C_90MEU5_P001_d$gene == 'CDKN2A',]

genelevel_s_C_K8F2FC_P001_d <- fread((here('data', 's_C_K8F2FC_P001_d__s_C_K8F2FC_N001_d.gene_level.txt')))
genelevel_s_C_K8F2FC_P001_d[genelevel_s_C_K8F2FC_P001_d$gene == 'CDKN2A',]

unique(cna_genelevel.txt$cn_state)

CDKN2A_cna <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'CDKN2A') %>% 
        select(sample, gene, cn_state)

CDKN2A_cna_all <- cna_genelevel.txt %>% 
        filter(gene == 'CDKN2A') 


CDKN2B_cna <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'CDKN2B') %>% 
        select(sample, gene, cn_state)  # s_C_5TP7FT_P001_d has CDKN2B CNLOH AFTER and CDKN2A DOUBLE LOSS AFTER

SETDB1_cna <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'SETDB1') %>% 
        select(sample, gene, cn_state)  

MLLT11_cna <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'MLLT11') %>% 
        select(sample, gene, cn_state)  ## exacly the same alterations as SETDB1 which makes sense

RB1_cna <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'RB1') %>% 
        select(sample, gene, cn_state) 

FRY_cna <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'FRY') %>% 
        select(sample, gene, chrom, gene_start, seg, segclust, gene_start, cn_state, filter)

FRY_cna_all <- cna_genelevel.txt %>% 
        group_by(sample) %>% 
        filter(gene == 'FRY') 

cna_chr13 <- cna_genelevel.txt %>% 
        filter(chrom == '13') 

cna_facets_run_info.txt$flags
###############################################################################################################

### Look for FRY mutations in CCLE to potentially validate these findings in vitro/silico
setwd("~/Dropbox (Partners HealthCare)/bootcamp/bootcamp_project/data/")

ccle_mut <- read.csv("CCLE_mutations.csv")
ccle_sample_info <- read.csv('sample_info.csv')

ccle_bladder_cells <- ccle_sample_info %>% filter(primary_disease == "Bladder Cancer") # 40 bladder cancer cell lines
FRY_ccle <- ccle_mut %>% filter(Hugo_Symbol == 'FRY') # 279 mutations

FRY_mut_celllines <- merge(FRY_ccle, ccle_sample_info, by = 'DepMap_ID')
FRY_urothelial <- FRY_mut_celllines %>% filter(primary_disease == "Bladder Cancer") # 3 cell lines, two primary samples (missense, silent), one met (Frame_Shift_Ins). UMUC13, UMUC7, UMUC9
ccle_urothelial <- ccle_mut %>% filter(primary_disease == "Bladder Cancer")
