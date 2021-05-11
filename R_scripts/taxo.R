source("script/r-stat_analysis/r_assign.r")

# Loading data
shared_path <- 
  "./data/N_pavm_mat_cluster_otu/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.RENAMED.shared.csv"
tax_path <-
  "./data/N_pavm_assign/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
shared <- assign_to_shared(shared_path = shared_path, tax_path = tax_path)
  

meta_path <- "./data/sample_metadata.tsv"
meta <- 
  read_tsv(file = meta_path, comment = "#")


# Little transformations
## R doesn't like colnames starting with a number
shared %>% 
  mutate(Group = str_replace(string = Group, 
                             pattern = "(.*)(.{1})$", 
                             replacement = "\\2\\1")) -> shared

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


######################################
taxo_level <- quo(Family)
taxo_level <- quo(!!as.name("Family"))
min_perc <- 5/100

shared %>%
  filter(otu_count > 0) %>% 
  select(Group, OTU, otu_count, group_count, !! (taxo_level)) %>% 
  group_by(!! (taxo_level), Group) %>%
  mutate(taxo_count = sum(otu_count)) %>%
  mutate(taxo_prop_group = taxo_count / group_count) %>%
  slice(1) %>% 
  filter(taxo_prop_group > min_perc) %>%
  select(-otu_count) %>% 
  ungroup() %>%
  group_by(Group) %>%
  arrange(desc(taxo_prop_group)) %>% 
  mutate(cum_prop = cumsum(taxo_prop_group)) %>%
  do(add_row(.data = .)) %>%
  ungroup() %>% 
  mutate_at(.vars = c("Group", "group_count"), 
            .funs = function(x) if_else(condition = is.na(x), 
                                                            true = lag(x), 
                                                            false = x)) %>%
  mutate_at(.vars = c("OTU", as_label(taxo_level)), 
            .funs = function(x) if_else(condition = is.na(x), 
                                       true = "Other", 
                                       false = x)) %>%
  group_by(Group) %>%
  mutate(taxo_count = if_else(condition = is.na(taxo_count), 
                             true = max(group_count) - sum(taxo_count, na.rm = TRUE), 
                             false = taxo_count)) %>%
  mutate(taxo_prop_group = taxo_count / group_count) %>%
  mutate(cum_prop = if_else(condition = is.na(cum_prop), 
                              true = lag(cum_prop) + taxo_prop_group, 
                              false = cum_prop)) %>%
  ungroup() %>% 
  left_join(y = meta, by = c("Group" = "sample")) %>%
  filter(ctrl != "ctrl_NEG") %>%
  group_by(patient) %>% 
  mutate(patient_num = group_indices()) %>%
  identity() -> shared_taxo_prop

#########################################
shared_taxo_prop %>% 
  filter(PAVM == FALSE, ctrl == "sample_POS") %>%
  arrange(patient) %>% 
  identity() %>% 
ggplot() +
    geom_bar(mapping = aes(x = reorder(Group, patient_num), 
                           y = taxo_prop_group, 
                           fill = !! (taxo_level),
                           group = Group), stat = "identity")

shared_taxo_prop %>% 
  filter(PAVM == TRUE, ctrl == "sample_POS") %>%
  arrange(patient) %>% 
  identity() %>% 
  ggplot() +
  geom_bar(mapping = aes(x = reorder(Group, patient_num), 
                         y = taxo_prop_group, 
                         fill = !! (taxo_level),
                         group = Group), stat = "identity")

shared_taxo_prop %>% 
  filter(ctrl == "sample_NEG") %>%
  arrange(patient) %>% 
  identity() %>% 
  ggplot() +
  geom_bar(mapping = aes(x = reorder(Group, patient_num), 
                         y = taxo_prop_group, 
                         fill = !! (taxo_level),
                         group = Group), stat = "identity")

shared_taxo_prop %>% 
  filter(ctrl == "ctrl_POS") %>%
  arrange(patient) %>% 
  identity() %>% 
  ggplot() +
  geom_bar(mapping = aes(x = reorder(Group, patient_num), 
                         y = taxo_prop_group, 
                         fill = !! (taxo_level),
                         group = Group), stat = "identity")


