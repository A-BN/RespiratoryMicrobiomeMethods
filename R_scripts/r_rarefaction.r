# Rarefaction curves PAVM, first batch
# antoine.bridier-nahmias@inserm.fr
# 2018-11-29

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


make_my_curve <-
  function(data, log = FALSE){
    namos <- deparse(substitute(data))
    data_curve <-
      data %>% 
      select(numsampled, contains("0.03")) %>% 
      gather(key = sample_name, value = count, -numsampled) %>% 
      mutate(log = log) %>% 
      mutate(numsampled = ifelse(test = log, yes = log10(numsampled), no = numsampled )) %>%
      ggplot() +
      geom_line(mapping = aes(x = numsampled, 
                              y = count, 
                              group = sample_name)) +
      if (log) {
        ggtitle( paste0("Log10_", namos))
      } else ggtitle(namos)
      
    return(data_curve)  
  }

plot_my_alpha <-
  function(data, rare_threshold, plot_title = ""){

    data_ <-
    data %>%
      mutate(closest = abs(numsampled - rare_threshold)) %>%
      filter(closest == min(closest)) 
    
    threshold_label <-
      format(x = data_$numsampled, scientific = TRUE)
    
    data_ <-
      data_ %>%
      select(-closest) %>% 
      select(-numsampled) %>%
      gather(key = sample_name, value = count) %>%
      separate(col = sample_name, into = c("measure", "sample"), sep = "-") %>% 
      mutate(measure = ifelse(test = measure == "0.03", yes = "mean", no = measure)) %>% 
      spread(key = measure, value = count) %>%
      replace(is.na(.), 0) %>% 
      identity()
    
    alpha_plot <-
    data_ %>% 
      ggplot(mapping = aes(x = sample)) + 
      geom_bar(mapping = aes(y = mean), stat = "identity", fill = "gray") +
      geom_errorbar(mapping = aes(ymin = lci, ymax = hci), alpha = 0.5) +
      ylab("") +
      ggtitle(plot_title) +
      coord_flip()
alpha_plot <- 
  ggdraw(
    add_sub(plot = alpha_plot, 
            label = paste("Rarefaction threshold =",
                          threshold_label)))
  
    return(list(data_,alpha_plot) )
  }

samples_left <- function(data){
  data_left <-
    data %>% 
    select(numsampled, contains("0.03")) %>% 
    gather(key = sample_name, value = count, -numsampled) %>% 
    group_by(sample_name) %>% 
    filter(!is.na(count)) %>% 
    mutate(max = max(numsampled)) %>% 
    slice(1) %>% 
    select(sample_name, max) %>% 
    ungroup() %>%
    arrange(desc(max)) %>% 
    identity()
  return(data_left)
}



##########################################
R_rarefa <- 
  read_tsv(file = "./data/R_pavm_alpha_1000/R_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.filter.groups.rarefaction")
R_shannon <- 
  read_tsv(file = "./data/R_pavm_alpha_1000/R_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.filter.groups.r_shannon")
# R_chao <- 
#   read_tsv(file = "./data/nofilter_R_pavm_alpha/R_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.groups.r_chao")
# R_coverage <- 
#   read_tsv(file = "./data/nofilter_R_pavm_alpha/R_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.groups.r_coverage")

make_my_curve(R_rarefa)
make_my_curve(R_shannon)
make_my_curve(R_shannon, log = TRUE)
# make_my_curve(R_chao)
# make_my_curve(R_coverage)

samples_left(R_rarefa) %>% View("R_sample_left")

R_threshold <- 5e3

R_notu <- 
  plot_my_alpha(data = R_rarefa, rare_threshold = R_threshold, plot_title = "R-primer, nOTU")[[2]]
R_Shannon <- 
  plot_my_alpha(data = R_shannon, rare_threshold = R_threshold, plot_title = "R-primer, Shannon's H'")[[2]]
cowplot::plot_grid(R_notu, R_Shannon)
plot(plot_my_alpha(data = R_shannon, rare_threshold = R_threshold)[[1]]$mean, plot_my_alpha(data = R_rarefa, rare_threshold = R_threshold)[[1]]$mean)

#################################################

N_rarefa <- 
  read_tsv(file = "./data/N_pavm_alpha_1000/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.filter.groups.rarefaction")
N_shannon <- 
  read_tsv(file = "./data/N_pavm_alpha_1000/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.0.03.filter.groups.r_shannon")
# N_chao <- 
#   read_tsv(file = "./data/nofilter_N_pavm_alpha/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.groups.r_chao")
# N_coverage <- 
#   read_tsv(file = "./data/nofilter_N_pavm_alpha/N_pavm_.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.groups.r_coverage")

make_my_curve(N_rarefa)
make_my_curve(N_shannon, log = TRUE)
# make_my_curve(N_chao)
# make_my_curve(N_coverage, TRUE)


samples_left(N_rarefa) %>% View("N_sample_left")

N_threshold <- 5e3

N_notu <- 
  plot_my_alpha(data = N_rarefa, rare_threshold = N_threshold, plot_title = "N-primer, nOTU")[[2]]
N_Shannon <- 
  plot_my_alpha(data = N_shannon, rare_threshold = N_threshold, plot_title = "N-primer, Shannon's H'")[[2]]
cowplot::plot_grid(N_notu, N_Shannon)

plot(plot_my_alpha(data = N_shannon, rare_threshold = N_threshold)[[1]]$mean, plot_my_alpha(data = N_rarefa, rare_threshold = N_threshold)[[1]]$mean)

cowplot::plot_grid(N_notu, R_notu, N_Shannon, R_Shannon)
