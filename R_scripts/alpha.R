# alpha diversity


library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)

# Functions
nobs <- function(x){
  sum(x > 0)
}

shannonH <- function(x, dig = 2){
  xln <- log(x)
  Hs <- -1 * sum(x * xln, na.rm = TRUE)
  Hs <- round(Hs, dig)
}

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



# Alpha table
alpha_df <-
  data.frame(sample = sample_names, 
             nobs = sapply(X = shared_prop, FUN = nobs),
             shannon = sapply(X = shared_prop, FUN = shannonH), 
             stringsAsFactors = FALSE)

alpha_df %>% 
  # mutate(patient = str_replace(string = sample, 
  #                              pattern = "[0-9]+([A-Z]+)N{1}", 
  #                              replacement = "\\1")) %>%
  left_join(y = meta, by = "sample") %>% 
  filter(sample != "N13BOM") %>% # Only LBA
  group_by(time, patient) %>% 
  mutate(ID = group_indices()) -> alpha_df

# Plotting time!
alpha_df2 <- 
  alpha_df %>%
  filter(startsWith(x = ctrl, "sample_")) %>%
  mutate(RT = as.factor(RT))#, PAVM = as.numeric(PAVM))

time.lab <- c("J-init", "J-VAP")
names(time.lab) <- c("J-3", "J6")
ctrl.lab <- c(sample_NEG = "control", sample_POS = "patient")

ggplot(data = alpha_df2, aes(group = ID)) +
  geom_line(mapping = aes(y = shannon, x = RT, color = patient)) +
  geom_point(mapping = aes(y = shannon,
                           x = RT,
                           color = patient),
             size = 5) +
  facet_wrap(~ time + ctrl, labeller = labeller(time = time.lab, ctrl = ctrl.lab))

ggplot(data = alpha_df2, aes(group = ID)) +
  geom_line(mapping = aes(y = shannon, x = time, color = patient)) +
  geom_point(mapping = aes(y = shannon,
                           x = time,
                           color = patient),
             size = 5) +
  facet_wrap(~ RT + ctrl)

ggplot(data = alpha_df2, aes(group = ID)) +
  geom_line(mapping = aes(y = nobs, x = RT, color = patient)) +
  geom_point(mapping = aes(y = nobs,
                           x = RT,
                           color = patient),
             size = 5) +
  facet_wrap(~ time + ctrl, labeller = labeller(time = time.lab, ctrl = ctrl.lab))


