# Centrifuge result for PAVM 2018-11 reads
# antoine.bridier-nahmias@inserm.fr
# 2018-12-05

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)

# Abundance in the reports is taking in account the genome size, we want to know
# the relative abundance in the sequences. The problem is that numReads is not
# clearly explained in Centrifuge's manual.

load_centrifuge <-
  function(report_path){
    report <-
      read_delim(report_path, 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
    total_unique <- sum(report$numUniqueReads)
    report <-
      report %>% 
        mutate(prop_unique = round(numUniqueReads / total_unique, digits = 3))  %>%
        identity()
    return(report)
  }

##### Read counts from X_pavm_.trim.contigs.good.unique.fasta

N_read_count <- as.numeric(
  str_replace(
    string = system(command = "wc -l ../../data/N_pavm_contigs/N_pavm_.trim.contigs.good.unique.fasta", 
                    intern = TRUE),
    pattern = "(^[0-9]+) .*",
    replacement = "\\1")) /2 # fasta format: each sequence has a header of 1 line 

R_read_count <- as.numeric(
  str_replace(
    string = system(command = "wc -l ../../data/R_pavm_contigs/R_pavm_.trim.contigs.good.unique.fasta", 
                    intern = TRUE),
    pattern = "(^[0-9]+) .*",
    replacement = "\\1")) /2 # fasta format: each sequence has a header of 1 line 

B_read_count <- as.numeric(
  str_replace(
    string = system(command = "wc -l ../../data/B_pavm_contigs/B_pavm_.trim.contigs.good.unique.fasta", 
                    intern = TRUE),
    pattern = "(^[0-9]+) .*",
    replacement = "\\1")) /2 # fasta format: each sequence has a header of 1 line 

#### Path to files
path_B_centrifuge <- 
  "~/Documents/budo/data/2018_11-PAVM/script/B_logs/centrifuge_report.tsv"

path_N_centrifuge <- 
  "~/Documents/budo/data/2018_11-PAVM/script/N_logs/centrifuge_report.tsv"

path_R_centrifuge <- 
  "~/Documents/budo/data/2018_11-PAVM/script/R_logs/centrifuge_report.tsv"


### Loading and transforming

B_report <-
  load_centrifuge(report_path = path_B_centrifuge) %>% 
  mutate(prop = round(numReads / B_read_count, digits = 3)) %>% 
  arrange(desc(prop)) %>% 
  identity()

N_report <-
  load_centrifuge(report_path = path_N_centrifuge)%>% 
  mutate(prop = round(numReads / N_read_count, digits = 3)) %>% 
  arrange(desc(prop)) %>% 
  identity()

R_report <-
  load_centrifuge(report_path = path_R_centrifuge)%>% 
  mutate(prop = round(numReads / R_read_count, digits = 3)) %>% 
  arrange(desc(prop)) %>% 
  identity()

write_delim(x = B_report, path = "../../data/reads/qc/B_metage.csv", delim = "\t", col_names = TRUE)
write_delim(x = N_report, path = "../../data/reads/qc/N_metage.csv", delim = "\t", col_names = TRUE)
write_delim(x = R_report, path = "../../data/reads/qc/R_metage.csv", delim = "\t", col_names = TRUE)

