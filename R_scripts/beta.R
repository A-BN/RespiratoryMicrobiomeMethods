# beta diversity

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)


# Loading data
shared_path <- 
  "./data/N_pavm_mat_cluster_otu/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.RENAMED.shared.csv"

shared <- 
  read_delim(shared_path, 
             "\t", escape_double = FALSE, trim_ws = TRUE)

meta_path <- "./data/sample_metadata.tsv"
meta <- 
  read_tsv(file = meta_path, comment = "#")


# Little transformations

## R doesn't like colnames starting with a number
sample_names <-   str_replace(string = shared$Group, 
                              pattern = "(.*)(.{1})$", 
                              replacement = "\\2\\1")

shared %>% 
  select(-label, -Group, -numOtus) %>% 
  t() %>% 
  as.data.frame() -> shared
colnames(shared) <- sample_names





meta %>% 
  mutate(RT = if_else(condition = (URT == 1), 
                      true = "upper", 
                      false = "lower")) %>% 
  mutate(time = if_else(condition = (time == 1), 
                        true = "J-3", 
                        false = "J6")) %>%
  # mutate(ctrl = if_else(condition = (ctrl == 0), 
  #                       true = "sample", 
  #                       false = FALSE)) %>%
  mutate(ctrl = case_when(
    ctrl == 0 ~ "sample_POS",
    ctrl == 1 ~ "sample_NEG",
    ctrl == 2 ~ "ctrl_NEG",
    ctrl == 3 ~ "ctrl_POS")) %>% 
  mutate(sample = str_replace(string = sample, 
                              pattern = "(.*)(.{1})$", 
                              replacement = "\\2\\1")) %>% 
  mutate(patient = str_replace(string = sample, 
                               pattern = "[NR]{1}[0-9]+([A-Z]+)", 
                               replacement = "\\1")) %>%
  # group_by(patient) %>% 
  # mutate(PAVM = if_else(condition = (sum(PAVM > 0) > 1),
  #                       true = TRUE,
  #                       false = FALSE)) %>%
  mutate(PAVM = if_else(condition = (PAVM == 1),
                        true = TRUE,
                        false = FALSE)) %>%
  # ungroup() %>% 
  select(-URT, -LRT) -> meta



# Shared file EDA
summary(rowSums(shared))
colSums(shared)

# Shared cleaning
n_reads_threashold <- 5
shared <- shared[rowSums(shared) > 5,]

# Shared proportions
shared_prop <- 
  as.data.frame(
    sapply(X = shared, FUN = function(x) x/sum(x)))

# Rarefied from mothur
# bray_df <- 
#   read_delim("data/N_pavm_beta_1000/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.braycurtis.0.03.column.ave.dist", 
#              "\t", escape_double = FALSE, col_names = c("sample_1", "sample_2", "dist"), 
#              trim_ws = TRUE)




























# beta div
beta_prop <- 
  t(shared_prop)

not_R <- stringr::str_starts(string = meta$sample, pattern = "R", negate = TRUE)
not_R_not_ctrl <- stringr::str_starts(string = meta$ctrl[not_R], pattern = "ctrl_", negate = TRUE)

beta_prop <-
  beta_prop[rownames(beta_prop) %in% meta$sample[not_R][not_R_not_ctrl], ]

# Bray curtis
beta_bray <- 
  vegan::vegdist(x = beta_prop, method = "bray", diag = TRUE)

bray_pcoa <- ape::pcoa(beta_bray)
biplot(bray_pcoa)
plot(x = bray_pcoa$vectors[,1], y = bray_pcoa$vectors[,2])

bray_df <-
  as.data.frame(as.matrix(beta_bray))

bray_df %>% 
  gather(key = "sample_1", value = "dist") %>% 
  mutate(sample_2 = rep(colnames(bray_df), ncol(bray_df))) %>% 
  left_join(y = meta, by = c("sample_1" = "sample")) %>% 
  mutate(PAVM_1 = PAVM, RT_1 = RT, ctrl_1 = ctrl, patient_1 = patient) %>% 
  mutate(PAVM = NULL, time = NULL, RT = NULL, patient = NULL, ctrl = NULL) %>% 
  left_join(y = meta, by = c("sample_2" = "sample")) %>% 
  mutate(PAVM_2 = PAVM, RT_2 = RT, ctrl_2 = ctrl, patient_2 = patient) %>% 
  mutate(PAVM = NULL, time = NULL, RT = NULL, patient = NULL, ctrl = NULL) -> bray_long

bray_long %>% 
  filter(sample_1 != sample_2) %>% 
  filter(ctrl_1 == ctrl_2) %>% 
  filter(PAVM_1 == PAVM_2) %>% 
  # half of the df is repeated
  rowwise() %>% 
  mutate(id = paste0(sort(c(sample_1, sample_2)), collapse = "")) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>%
  select(-id) %>%
  mutate(status = paste(ctrl_1, PAVM_1, sep = "_")) %>%
  mutate(status = case_when(status == "sample_NEG_FALSE" ~ "NEG",
                            status == "sample_POS_FALSE" ~ "PRE",
                            status == "sample_POS_TRUE" ~ "PAVM")) %>%
  mutate(status = factor(status, levels = c("NEG", "PRE", "PAVM"))) %>% 
  identity() -> bray_long_filt

ggplot(data = bray_long_filt) +
  geom_boxplot(mapping = aes(x = status, y = dist)) +
  facet_wrap(~RT_1)


# Unifrac
otu_tree <-
  ape::read.tree(file = "./data/N_pavm_tree/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.rep.tre")

GUniFrac::GUniFrac(otu.tab = t(shared), tree = otu_tree)




bray_disper <- vegan::betadisper(d = beta_bray, group = meta$ctrl[not_R][not_R_not_ctrl])
plot(bray_disper)
vegan::permutest(x = bray_disper, pairwise = TRUE)

vegan::adonis(formula = beta_bray ~ time, data = meta[not_R, ][not_R_not_ctrl, ])
