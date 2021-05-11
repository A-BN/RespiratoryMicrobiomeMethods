# assign taxonomy to shared file PAVM, first batch
# antoine.bridier-nahmias@inserm.fr
# 2018-12-05

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)

assign_to_shared <- 
  function(shared_path, tax_path){
    
    shared <- read_delim(shared_path, 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
    tax <- read_delim(tax_path, 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
    
    shared <-
      shared %>% 
      select(-label, -numOtus) %>%
      gather(key = OTU, value = count, -Group) %>% 
      identity()
    
    tax <- 
      tax %>%
      select(-Size) %>% 
      mutate(Taxonomy = str_replace_all(pattern = "\\([0-9]+\\)", replacement = "", string = Taxonomy)) %>% 
      # separate(col = Taxonomy, sep = ";", into = paste0("lvl", 1:6)) %>%
      separate(col = Taxonomy, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>% 
      identity()
    
    cat("Number of 'unclassified' by taxonomic rank: \n")
    tax %>% 
      select(-OTU) %>% 
      sapply(FUN = function(x) sum(str_detect(string = x, pattern = "unclassified"))) %>% 
      cat("\n")
    
    shared <- 
      left_join(x = shared, y = tax, by = "OTU") %>%
      rename(otu_count = count) %>% 
      group_by(Group) %>% 
      mutate(group_count = sum(otu_count)) %>% 
      ungroup() %>% 
      identity()
    
    return(shared)
  }

largest_otu_plot <-
  function(shared_tax = N_shared, n_otu = 10, tax_level = "Genus"){
    
    shared_plot <-
      shared_tax %>% 
      group_by(Group) %>%
      arrange(desc(otu_count)) %>% 
      slice(1:n_otu) %>% 
      ungroup() %>%
      mutate(otu_prop = otu_count / group_count) %>% 
      identity()
    
    plotos <-
      shared_plot %>%
      mutate(taxo = !!ensym(tax_level)) %>% 
      group_by(Group) %>%
      do( plots =
            ggplot(data = ., mapping = aes(x = taxo, y = otu_prop, group = .$taxo)) +
            geom_bar(stat = "identity", fill = "gray") +
            ggtitle(.$Group) +
            # xlab("") +
            coord_flip()
          )
    return(list(plotos, shared_plot))
  }


# shared_path <-
#   "~/Documents/budo/data/2018_11-PAVM/data/N_pavm_mat_cluster_otu/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.shared"
# 
# tax_path <-
#   "~/Documents/budo/data/2018_11-PAVM/data/N_pavm_assign/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"



# 
# shared_N_path <- 
#   "./data/N_pavm_mat_cluster_otu/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.shared"
# 
# tax_N_path <- 
#   "./data/N_pavm_assign/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
# 
# N_shared <-
#   assign_to_shared(shared_path = shared_N_path, 
#                    tax_path = tax_N_path)
# N_largest <-
#   largest_otu_plot(shared_tax = N_shared, n_otu = "50", tax_level = "Genus")[[1]]
# ############"
# shared_R_path <- 
#   "./data/R_pavm_mat_cluster_otu/R_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.shared"
# 
# tax_R_path <- 
#   "./data/R_pavm_assign/R_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
# 
# R_shared <-
#   assign_to_shared(shared_path = shared_R_path, 
#                    tax_path = tax_R_path)
# R_largest <-
#   largest_otu_plot(shared_tax = R_shared, n_otu = 100, tax_level = "Genus")[[1]]
# 
# ################
# 
# # Merging of the "X_largest dataframe"
# 
# merge_largest <-
# N_largest %>% 
#   # rbind(R_largest) %>% 
#   mutate(patient = str_replace(string = Group, pattern = "(.*).{1}$", replacement = "\\1"))
# merge_largest$index <- 1:nrow(merge_largest)
# 
# 
# for (i in 1:length(unique(merge_largest$patient))) {
#   # i <- 13
#   curr_patient <-
#     (merge_largest$patient)[i]
#   curr_index <- 
#     merge_largest$patient == curr_patient
#   curr_plot <- 
#     plot_grid(plotlist = merge_largest$plots[curr_index])
#   
#   # ggsave(plot = curr_plot, 
#   #        filename = paste0(i, ".pdf"), 
#   #        path = "../../2019-01-17_for_MelFro/", 
#   #        width = 40, 
#   #        height = 5, 
#   #        units = "cm")
#  }
# 
# 
# sapply(X = N_largest$plots, FUN = function(x) plot(x) )
# 
# ####### TODO
# # Add a color option to largest_otu_plot()

