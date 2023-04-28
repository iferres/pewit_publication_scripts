setwd("/mnt/3tb/iferres/pewit_publication_scripts/")

library(magrittr)
library(future)
library(future.apply)
plan(multisession, workers = 6, gc = TRUE)

###################
# 0. Dependencies #
###################
#' Install R packages to test, and download containers with external software to
#' test.

##  0.1 R-packages

# 0.1.1 PEWIT
if (!require(devtools)) install.packages("devtools")
devtools::install_github('iferres/pewit', ref = '97da040')


# 0.1.2 MICROPAN
install.packages('micropan')
packageVersion("micropan")
# [1] ‘2.1’

# 0.2 Setup containers
containers <- future(source("Helper_scripts/0.2_Download-containers.R", local = TRUE)$value)


##################
# 1. Simulations #
##################

## 1.1 Simulate pangenomes

dir.create("1.1_Run-simulations")
source("Helper_scripts/1.1_Run-simulations.R")
nes <-  c(2e+09, 1e+10, 5e+10, 25e+10, 125e+10)
rep <- 5L
set.seed(123)
simulations <- expand.grid(nes, seq_len(rep)) %>%
  setNames(c("ne", "rep")) %>% {
    future_Map(function(n, r){
      
      simulate_pangenome(ne = n, repli = r)
      
    }, n = .$ne, r = .$rep, 
    future.stdout = FALSE,
    future.seed = TRUE)
  }


## 1.2 Format simulations

sim_fastas <- list.files(path = "1.1_Run-simulations/", 
                         pattern = "genome\\d+[.]fasta$", 
                         recursive = TRUE, 
                         full.names = TRUE)


# 1.2.1 Translate

sim_faa <- future_lapply(sim_fastas, function(x){
  out <- sub('[.]fasta$', '.faa', x)
  ffn <- read.fasta(x)
  faas <- lapply(ffn, translate)
  write.fasta(faas, names = names(faas), file.out = out)
  return(out)
}, future.packages = "seqinr")

source("Helper_scripts/1.2_Format-simulations.R")



# 1.2.2 Generate gffs

sim_gff <- future_lapply(sim_fastas, ffn2gff, future.packages = "seqinr")


# 1.2.3 Generate embl using gff3toembl

sim_embl <- future_lapply(sim_gff, gff3_2_embl, container = value(containers)$gff3toembl, future.stdout = FALSE)


# 1.2.4 Generate gbks

sim_gbk <- future_lapply(sim_embl, embl_2_gbk, future.stdout = FALSE)



## 1.3 Run clustering benchmarks
sim_dirs <- dir(path = "1.1_Run-simulations/", full.names = T)

## 1.3.1 ROARY i70

source("Helper_scripts/1.3_Run-Clustering-Benchmarks.R")
plan(multisession, workers = 6, gc = TRUE)
roary_i70 <- future_lapply(sim_dirs, function(x){
  bench_roary(
    din = x,
    dout = 'roaryi70_vs_sim', 
    name = "roary_i70",
    container = value(containers)$roary,
    container_exec_args = "-B $PWD",
    cmd = "roary",
    cmd_args = "-p 1 -i 70",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE, 
future.conditions = character(0))
future:::ClusterRegistry("stop")

roary_i70_df <- roary_i70 %>%
  do.call(rbind, .)
write.table(roary_i70_df, file = "1.3_Run-clustering-benchmark/roaryi70_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)
  

## 1.3.2 ROARY i95

plan(multisession, workers = 6, gc = TRUE)
roary_i95 <- future_lapply(sim_dirs, function(x){
  bench_roary(
    din = x,
    dout = 'roaryi95_vs_sim', 
    name = "roary_i95",
    container = value(containers)$roary,
    container_exec_args = "-B $PWD",
    cmd = "roary",
    cmd_args = "-p 1 -i 95",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE, 
future.conditions = character(0))
future:::ClusterRegistry("stop")

roary_i95_df <- roary_i95 %>%
  do.call(rbind, .)
write.table(roary_i95_df, file = "1.3_Run-clustering-benchmark/roaryi95_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)

# 1.3.3 PANX

plan(multisession, workers = 6, gc = TRUE)
panx <- future_lapply(sim_dirs, function(x){
  bench_panx(
    din = x, 
    dout = 'panx_vs_sim', 
    name = "panx",
    container = value(containers)$panx,
    container_exec_args = "-B $PWD",
    cmd = "panX.py",
    cmd_args = "-sl Simulated_bacteria -cg 0.7 -st 1 2 3 4 5 6 -t 1 -ct",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE,
future.conditions = character(0))
future:::ClusterRegistry("stop")

panx_df <- panx %>%
  do.call(rbind, .)
write.table(panx_df, file = "1.3_Run-clustering-benchmark/panx_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)

# 1.3.4 MICROPAN COMPLETE 

plan(multisession, workers = 6, gc = TRUE)
micropanBlast_Complete <- future_mapply(function(x, job){
  bench_micropan_blast(din = x, 
                       dout = "micropanBlastComplete_vs_sim", 
                       name = "micropanBlast_Complete", 
                       linkage = "complete", 
                       job = job,
                       work_dir = "1.3_Run-clustering-benchmark")},
  x = sim_dirs, job = seq_along(sim_dirs),
  future.stdout = FALSE,
  future.conditions = character(0), 
  future.packages = "micropan")
future:::ClusterRegistry("stop")

micropanBlast_Complete_df <- micropanBlast_Complete %>%
  t()
write.table(micropanBlast_Complete_df, file = "1.3_Run-clustering-benchmark/micropanBlastComplete_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)

# 1.3.5 MICROPAN AVERAGE

plan(multisession, workers = 6, gc = TRUE)
micropanBlast_Average <- future_mapply(function(x, job){
  bench_micropan_blast(din = x, 
                       dout = "micropanBlastAverage_vs_sim", 
                       name = "micropanBlast_Average", 
                       linkage = "average", 
                       job = job,
                       work_dir = "1.3_Run-clustering-benchmark")},
  x = sim_dirs, job = seq_along(sim_dirs),
  future.stdout = FALSE,
  future.conditions = character(0), 
  future.packages = "micropan")
future:::ClusterRegistry("stop")

micropanBlast_Average_df <- micropanBlast_Average %>%
  t()
write.table(micropanBlast_Average_df, file = "1.3_Run-clustering-benchmark/micropanBlastAverage_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)

# 1.3.6 MICROPAN SINGLE

plan(multisession, workers = 6, gc = TRUE)
micropanBlast_Single <- future_mapply(function(x, job){
  bench_micropan_blast(din = x, 
                       dout = "micropanBlastSingle_vs_sim", 
                       name = "micropanBlast_Single", 
                       linkage = "single", 
                       job = job,
                       work_dir = "1.3_Run-clustering-benchmark")},
  x = sim_dirs, job = seq_along(sim_dirs),
  future.stdout = FALSE,
  future.conditions = character(0), 
  future.packages = "micropan")
future:::ClusterRegistry("stop")

micropanBlast_Single_df <- micropanBlast_Single %>% t()
write.table(micropanBlast_Single_df, file = "1.3_Run-clustering-benchmark/micropanBlastSingle_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)

# 1.3.7 PANAROO STRICT

plan(multisession, workers = 6, gc = TRUE)
panaroo_strict <- future_lapply(sim_dirs, function(x){
  bench_panaroo(
    din = x, 
    dout = "panarooStrict_vs_sim", 
    name = "panarooStrict",
    container = value(containers)$panaroo,
    container_exec_args = "-B $PWD",
    cmd = "panaroo",
    cmd_args = "--clean-mode strict --threads 1",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE)
future:::ClusterRegistry("stop")

panaroo_strict_df <- panaroo_strict %>%
  do.call(rbind, .)
write.table(panaroo_strict_df, file = "1.3_Run-clustering-benchmark/panarooStrict_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)


# 1.3.8 PANAROO MODERATE

plan(multisession, workers = 6, gc = TRUE)
panaroo_moderate <- future_lapply(sim_dirs, function(x){
  bench_panaroo(
    din = x, 
    dout = "panarooModerate_vs_sim", 
    name = "panarooModerate",
    container = value(containers)$panaroo,
    container_exec_args = "-B $PWD",
    cmd = "panaroo",
    cmd_args = "--clean-mode moderate --threads 1",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE)
future:::ClusterRegistry("stop")

panaroo_moderate_df <- panaroo_moderate %>%
  do.call(rbind, .)
write.table(panaroo_moderate_df, file = "1.3_Run-clustering-benchmark/panarooModerate_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)


# 1.3.9 PANAROO SENSITIVE

plan(multisession, workers = 6, gc = TRUE)
panaroo_sensitive <- future_lapply(sim_dirs, function(x){
  bench_panaroo(
    din = x, 
    dout = "panarooSensitive_vs_sim", 
    name = "panarooSensitive",
    container = value(containers)$panaroo,
    container_exec_args = "-B $PWD",
    cmd = "panaroo",
    cmd_args = "--clean-mode sensitive --threads 1",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE)
future:::ClusterRegistry("stop")

panaroo_sensitive_df <- panaroo_sensitive %>%
  do.call(rbind, .)
write.table(panaroo_sensitive_df, file = "1.3_Run-clustering-benchmark/panarooSensitive_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)


# 1.3.10 PEWIT

plan(multisession, workers = 6, gc = TRUE)
pewit <- future_lapply(sim_dirs, function(x){
  bench_pewit(
    din = x, 
    dout = "pewit_vs_sim", 
    name = "pewit",
    hmm = "Data/Pfam/Pfam-A.hmm",
    dat = "Data/Pfam/Pfam-A.hmm.dat",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE)
future:::ClusterRegistry("stop")

pewit_df <- pewit %>%
  do.call(rbind, .)
write.table(pewit_df, file = "1.3_Run-clustering-benchmark/pewit_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)


# 1.3.11 PEWIT No Pfam

plan(multisession, workers = 6, gc = TRUE)
pewit_noPfam <- future_lapply(sim_dirs, function(x){
  bench_pewit(
    din = x, 
    dout = "pewitNoPfam_vs_sim", 
    name = "pewitNoPfam",
    work_dir = "1.3_Run-clustering-benchmark"
  )
}, 
future.stdout = FALSE)
future:::ClusterRegistry("stop")

pewit_noPfam_df <- pewit_noPfam %>%
  do.call(rbind, .)
write.table(pewit_noPfam_df, file = "1.3_Run-clustering-benchmark/pewitNoPfam_vs_sim_results.tsv", quote = FALSE, sep = "\t", row.names = F)




## 1.4 Analysis

dir.create("1.4_Analysis_Clustering")

sim_results <- list.files(path = "1.3_Run-clustering-benchmark/", 
                          pattern = "_vs_sim_results.tsv", 
                          full.names = T)

names(sim_results) <- sim_results %>%
  sub(".+//", "", .) %>%
  sub("_vs_sim.+", "", .)

sim_df <- sim_results %>%
  lapply(read.csv, sep="\t") %>%
  do.call(rbind, .)

software_col <- c("pewit" = "#DC0000FF", 
                   "pewitNoPfam" = "#E64B35FF",
                   "roary_i70" = "#8F7700FF",
                   "roary_i95" = "#DF8F44FF",
                   "panx" = "#00A087FF",
                   "micropanBlast_Complete" = "#003C67FF",
                   "micropanBlast_Average" = "#3C5488FF",
                   "micropanBlast_Single" = "#8491B4FF",
                   "panarooStrict" = "#CD534CFF",
                   "panarooModerate" = "#7E6148FF",
                   "panarooSensitive" = "#B09C85FF")

swlabels <- c(pewit = "Pewit", 
              pewitNoPfam = "Pewit\nwithout Pfam", 
              roary_i70 = "Roary -i 70", 
              roary_i95 = "Roary -i 95", 
              panx = "PanX\nDIAMOND All vs All", 
              micropanBlast_Complete = "Micropan\nBlastp All vs All\n Complete", 
              micropanBlast_Average = "Micropan\nBlastp All vs All\n Average", 
              micropanBlast_Single = "Micropan\nBlastp All vs All\n Single", 
              panarooStrict = "Panaroo\nClean mode:\nStrict", 
              panarooModerate = "Panaroo\nClean mode:\nModerate", 
              panarooSensitive = "Panaroo\nClean mode:\nSensitive")


##################################################
## Compute all binary classification statistics ##
##################################################
sim_df %<>%
  within(Recall <- TP / (TP + FN)) %>% # aka Sensitivity
  within(Precision <- TP / (TP + FP)) %>%
  within(Specificity <- TN / (TN + FP)) %>%
  within(Accuracy <- (TP + TN) / (TP + TN + FP + FN)) %>%
  
  within(TPR <- TP / (TP + FN)) %>% # True Positive Rate
  within(TNR <- TN / (TN + FP)) %>% # True Negative Rate
  within(PPV <- TP / (TP + FP)) %>% # Positive Prediction Value
  within(NPV <- TN / (TN + FN)) %>% # Negative Prediction Value
  
  within(FPR <- FP / (FP + TN)) %>% # False Positive Rate
  within(FNR <- FN / (FN + TP)) %>% # False Negative Rate
  within(FOR <- FN / (FN + TN)) %>% # False Omission Rate
  within(FDR <- FP / (FP + TP)) %>% # False Discovery Rate
  
  within(F1_Score <- 2 * ((Precision * Recall) /
                            Precision + Recall)) %>%
  within(Youden_J <- Recall + Specificity - 1) %>%
  within(Fowlkes_Mallows <- sqrt( (TP / (TP + FP)) * (TP / (TP + FN)) ) ) %>%
  # Matthews Correlation Coeficient
  within(MCC <- sqrt(PPV * TPR * TNR * NPV) - 
           sqrt(FDR * FNR * FPR * FOR)) %>%
  within(Software <- factor(Software, levels = names(software_col)))

write.table(sim_df, file = "1.4_Analysis_Clustering/simulation_results.tsv",quote = FALSE, sep = "\t", row.names = FALSE)


#####################
## Plot Statistics ##
#####################
library(ggplot2)

# Recall
sim_df %>%
  ggplot(aes(y = Recall, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("Recall = ", frac("TP", "TP + FN"))))
ggsave("1.4_Analysis_Clustering/Recall.png")

# Precision
sim_df %>%
  ggplot(aes(y = Precision, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("Precision = ", frac("TP", "TP + FP"))))
ggsave("1.4_Analysis_Clustering/Precision.png")


# Specificity (True negative rate)
sim_df %>%
  ggplot(aes(y = Specificity, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("Specificity = ", frac("TN", "FP + TN"))))
ggsave("1.4_Analysis_Clustering/Specificity.png")

# Accuracy
sim_df %>%
  ggplot(aes(y = Accuracy, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("Accuracy = ", frac("TP + TN", "TP + FP + FP + FN"))))
ggsave("1.4_Analysis_Clustering/Accuracy.png")


# False Positive Rate
sim_df %>%
  ggplot(aes(y = FPR, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("False Positive Rate = ", frac("FP", "FP + TN"))))
ggsave("1.4_Analysis_Clustering/FPR.png")

# False Discovery Rate
sim_df %>%
  ggplot(aes(y = FDR, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("False Discovery Rate = ", frac("FP", "TP + FP"))))
ggsave("1.4_Analysis_Clustering/FDR.png")

# False Omission Rate
sim_df %>%
  ggplot(aes(y = FOR, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  ggtitle(expression(paste("False Omission Rate = ", frac("FN", "FN + TN"))))
ggsave("1.4_Analysis_Clustering/FOR.png")

# F1 Score
sim_df %>%
  ggplot(aes(y = F1_Score, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) + 
  
  ggtitle(expression(paste(F[1], "score = 2 .", frac("Precision . Recall", "Precision + Recall"))))
ggsave("1.4_Analysis_Clustering/F1_Score.png")

# Fowlkes-Mallows index
sim_df %>%
  ggplot(aes(y = Fowlkes_Mallows, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=6),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
    )
ggsave("1.4_Analysis_Clustering/Fowlkes_Mallows.png")


# Matthews Correlation Coeficient
f1a <- sim_df %>%
  ggplot(aes(y = MCC, color = Software, x = Software)) +
  geom_boxplot() + 
  scale_color_manual(values = software_col) +
  facet_wrap(~Ne, nrow = 1) + 
  theme_bw() + 
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    text = element_text(family = "serif")
  ) + 
  guides(color=FALSE);f1a
ggsave("1.4_Analysis_Clustering/Matthews-Correlation-Coefficient.png", height = 6)


# Precision vs Recall
f1b <- sim_df %>%
  ggplot(aes(x = Precision, y = Recall, color = Software)) +
  geom_point(aes(shape=as.factor(Ne)), size=5, alpha=.5) + 
  geom_abline(slope = 1, intercept = 0, linetype=2, size = .2) +
  facet_wrap(~Software, labeller = as_labeller(swlabels)) +
  scale_color_manual(values = software_col) +
  scale_fill_manual(values = software_col) +
  scale_shape_discrete(name = "Number of\ngenerations") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = .5, vjust = .5), 
        strip.background = element_rect(fill = NA), 
        text = element_text(family = "serif"));f1b
ggsave("1.4_Analysis_Clustering/Recall_vs_Precision.png", height = 6)
ggsave("1.4_Analysis_Clustering/Recall_vs_Precision.pdf", height = 6)

library(patchwork)
fig1 <- (f1a / f1b) + 
  patchwork::plot_layout(guides = "collect", heights = c(1,2)) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave("1.4_Analysis_Clustering/Figure2.png", plot = fig1, width = 13.3, height = 6)
ggsave("1.4_Analysis_Clustering/Figure2.pdf", plot = fig1, width = 25, height = 20, units = "cm")

sim_df %>%
  ggplot(aes(x = Precision, y = Recall, color = Software)) +
  geom_point(size=4, alpha=.5) + 
  geom_abline(slope = 1, intercept = 0, linetype=2, size = .2) +
  facet_wrap(~Ne, nrow = 1) +
  scale_color_manual(values = software_col) +
  scale_fill_manual(values = software_col) +
  # scale_shape_discrete(name = "Effective population\nsize") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = .5, vjust = .5), 
        strip.background = element_rect(fill = NA))
ggsave("1.4_Analysis_Clustering/Recall_vs_Precision_facet-Ne.png", width = 10, height = 4)




####################
# 2. Running times #
####################

# 2.1 List Datasets
cfetus_gffs <- list.files(path = "Data/Dataset_Cfetus/", 
                          pattern = "[.]gff$", 
                          full.names = TRUE)

cfetus_gbks <- list.files(path = "Data/Dataset_Cfetus/", 
                          pattern = "[.]gbk$", 
                          full.names = TRUE)

cfetus_faas <- list.files(path = "Data/Dataset_Cfetus/", 
                          pattern = "[.]gbk$", 
                          full.names = TRUE)

cfetus_embl <- list.files(path = "Data/Dataset_Cfetus/", 
                          pattern = "[.]embl$", 
                          full.names = TRUE)



# 2.2 Sample Datasets

NGen <- c(10, 20, 30, 50, 70, 90, 120, 150) # Number of genomes
reps <- 1:5                                 # Number of replicates

set.seed(123)
df <- lapply(NGen, function(x){
  lp <- lapply(reps, function(y){
    gffs <- sample(cfetus_gffs, x, replace = FALSE)
    cbind(gffs, x, y)
  })
  do.call(rbind, lp)
})
df <- as.data.frame(do.call(rbind, df), stringsAsFactors = FALSE)
colnames(df) <- c('gff_name', 'Num_genomes', 'replicate')
df$Num_genomes <- as.integer(df$Num_genomes)
df$replicate <- as.integer(df$replicate)
df$set_id <- paste0('N', df$Num_genomes, 'R', df$replicate)

write.table(df, file = 'Data/NumGenomes_set.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
# df <- read.csv("Data/NumGenomes_set.tsv", sep = "\t")


# 2.3 Running time benchmark
source("Helper_scripts/2.3_Runtime-Benchmarks.R")

# 2.3.1 Pewit

plan(multisession, workers = 3, gc = TRUE)
rt_pewit <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_pewit(
      gffs = x$gff_name,
      hmm_pfam = "Data/Pfam/Pfam-A.hmm",
      dat_pfam = "Data/Pfam/Pfam-A.hmm.dat"
    )
    data.frame(
      Software = "pewit",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_pewit, 
            file = "2.3_Runtime-Benchmarks/runtime_pewit.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.2 Pewit without Pfam

plan(multisession, workers = 3, gc = TRUE)
rt_pewitNoPfam <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_pewit(
      gffs = x$gff_name
    )
    data.frame(
      Software = "pewitNoPfam",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_pewitNoPfam, 
            file = "2.3_Runtime-Benchmarks/runtime_pewitNoPfam.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.3 Roary -i 70

plan(multisession, workers = 3, gc = TRUE)
rt_roary_i70 <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_roary(
      gffs = x$gff_name,
      set_id = x$set_id[1],
      cmd_args = "-p 1 -i 70",
      dout = "roary_i70"
    )
    data.frame(
      Software = "roary_i70",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_roary_i70, 
            file = "2.3_Runtime-Benchmarks/runtime_roary_i70.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.4 Roary -i 95

plan(multisession, workers = 3, gc = TRUE)
rt_roary_i95 <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_roary(
      gffs = x$gff_name,
      set_id = x$set_id[1],
      cmd_args = "-p 1 -i 95",
      dout = "roary_i95"
    )
    data.frame(
      Software = "roary_i95",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_roary_i95, 
            file = "2.3_Runtime-Benchmarks/runtime_roary_i95.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.5 PanX

plan(multisession, workers = 3, gc = TRUE)
rt_panx <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_panx(
      gbks = sub("gff$", "gbk", x$gff_name),
      set_id = x$set_id[1],
      dout = "panx"
    )
    data.frame(
      Software = "panx",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE, 
  future.conditions = character(0)) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_panx, 
            file = "2.3_Runtime-Benchmarks/runtime_panx.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.6 Panaroo strict mode

plan(multisession, workers = 3, gc = TRUE)
rt_panarooStrict <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_panaroo(
      gffs = x$gff_name,
      set_id = x$set_id[1],
      dout = "panarooStrict",
      cmd_args = "--clean-mode strict --threads 1"
    )
    data.frame(
      Software = "panarooStrict",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE, 
  future.conditions = character(0)) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_panarooStrict, 
            file = "2.3_Runtime-Benchmarks/runtime_panarooStrict.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.7 Panaroo moderate mode

plan(multisession, workers = 3, gc = TRUE)
rt_panarooModerate <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_panaroo(
      gffs = x$gff_name,
      set_id = x$set_id[1],
      dout = "panarooModerate",
      cmd_args = "--clean-mode moderate --threads 1"
    )
    data.frame(
      Software = "panarooModerate",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE, 
  future.conditions = character(0)) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_panarooModerate, 
            file = "2.3_Runtime-Benchmarks/runtime_panarooModerate.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.8 Panaroo sensitive mode

plan(multisession, workers = 3, gc = TRUE)
rt_panarooSensitive <- split(df, df$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_panaroo(
      gffs = x$gff_name,
      set_id = x$set_id[1],
      dout = "panarooSensitive",
      cmd_args = "--clean-mode sensitive --threads 1"
    )
    data.frame(
      Software = "panarooSensitive",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE, 
  future.conditions = character(0)) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_panarooSensitive, 
            file = "2.3_Runtime-Benchmarks/runtime_panarooSensitive.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.9 Micropan Blastp All vs All

plan(multisession, workers = 3, gc = TRUE)
rt_micropan_BlastAllAll <- df[which(df$Num_genomes<=50), ] %>%
  split(., .$set_id) %>%
  future_lapply(function(x){
    rt <- runtime_Micropan_BlastAllAll(
      faas = sub("gff$", "faa", x$gff_name),
      set_id = x$set_id[1],
      dout = "micropan_BlastAllAll"
    )
    data.frame(
      Software = "micropanBlastAllAll",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  },future.stdout = FALSE, 
  future.conditions = character(0)) %>%
  do.call(rbind, .) 
future:::ClusterRegistry("stop")
write.table(rt_micropan_BlastAllAll, 
            file = "2.3_Runtime-Benchmarks/runtime_micropanBlastAllAll.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.9.1 Micropan Clustering single linkage

plan(multisession, workers = 3, gc = TRUE)
rt_micropanBlast_Single <- df[which(df$Num_genomes<=50), ] %>%
  split(., .$set_id) %>% 
  future_lapply(function(x){
    rt <- runtime_Micropan_Clust(
      set_id = x$set_id[1],
      linkage = "single"
    )
    data.frame(
      Software = "micropanBlast_Single",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  }, future.stdout = FALSE,
  future.conditions = character(0)) %>%
  do.call(rbind, .)
future:::ClusterRegistry("stop")
write.table(rt_micropanBlast_Single, 
            file = "2.3_Runtime-Benchmarks/runtime_micropanBlastAllAll_Single.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.9.2 Micropan Clustering average linkage

plan(multisession, workers = 3, gc = TRUE)
rt_micropanBlast_Average <- df[which(df$Num_genomes<=50), ] %>%
  split(., .$set_id) %>% 
  future_lapply(function(x){
    rt <- runtime_Micropan_Clust(
      set_id = x$set_id[1],
      linkage = "average"
    )
    data.frame(
      Software = "micropanBlast_Average",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  }, future.stdout = FALSE,
  future.conditions = character(0)) %>%
  do.call(rbind, .)
future:::ClusterRegistry("stop")
write.table(rt_micropanBlast_Average, 
            file = "2.3_Runtime-Benchmarks/runtime_micropanBlastAllAll_Average.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)

# 2.3.9.3 Micropan Clustering complete linkage

plan(multisession, workers = 3, gc = TRUE)
rt_micropanBlast_Complete <- df[which(df$Num_genomes<=50), ] %>%
  split(., .$set_id) %>% 
  future_lapply(function(x){
    rt <- runtime_Micropan_Clust(
      set_id = x$set_id[1],
      linkage = "complete"
    )
    data.frame(
      Software = "micropanBlast_Complete",
      Num_genomes = x$Num_genomes[1],
      Replicate = x$replicate[1],
      set_id = x$set_id[1],
      elapsed_time = rt$system_time[["elapsed"]],
      exit_code = rt$exit
    )
  }, future.stdout = FALSE,
  future.conditions = character(0)) %>%
  do.call(rbind, .)
future:::ClusterRegistry("stop")
write.table(rt_micropanBlast_Complete, 
            file = "2.3_Runtime-Benchmarks/runtime_micropanBlastAllAll_Complete.tsv",
            quote = F, 
            sep = "\t", 
            row.names = FALSE)


# 2.4 Running-time analysis

rt_tsvs <- list.files(path = "2.3_Runtime-Benchmarks/", 
                      pattern = "[.]tsv",
                      full.names = TRUE) 
rt_tsvs <- rt_tsvs[-grep("BlastAllAll[.]tsv", rt_tsvs)]
dir.create("2.4_Analysis_Running-time")

library(ggplot2)

# list(rt_pewit,
#      rt_pewitNoPfam,
#      rt_roary_i70,
#      rt_roary_i95,
#      rt_panx,
#      rt_panarooStrict,
#      rt_panarooModerate,
#      rt_panarooSensitive,
#      rt_micropanBlast_Single,
#      rt_micropanBlast_Average,
#      rt_micropanBlast_Complete
#      ) %>%

swlabels <- c(pewit = "Pewit", 
              pewitNoPfam = "Pewit\nwithout Pfam", 
              roary_i70 = "Roary -i 70", 
              roary_i95 = "Roary -i 95", 
              panx = "PanX", 
              micropanBlast_Complete = "Micropan\nBlastp All vs All\n Complete", 
              micropanBlast_Average = "Micropan\nBlastp All vs All\n Average", 
              micropanBlast_Single = "Micropan\nBlastp All vs All\n Single", 
              panarooStrict = "Panaroo\nClean mode: Strict", 
              panarooModerate = "Panaroo\nClean mode: Moderate", 
              panarooSensitive = "Panaroo\nClean mode: Sensitive"
)

f2a <- rt_tsvs %>%
  lapply(read.csv, sep="\t", header = TRUE) %>%
  do.call(rbind, .) %>%
  subset(Num_genomes < 150) %>%
  within(Software <- factor(Software, levels = names(software_col))) %>%
  ggplot(aes(x=Num_genomes, y=elapsed_time)) + 
  geom_point(aes(fill=Software), shape = 21, alpha = .3) + 
  stat_summary(fun = mean, geom = "line", mapping = aes(color=Software, group=Software)) +
  scale_color_manual(values = software_col) +
  scale_fill_manual(values = software_col) +
  scale_x_continuous(breaks = NGen) +
  # gghighlight::gghighlight() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = .5, vjust = .5), 
        strip.background = element_rect(fill = NA),
        text = element_text(family = "serif")) + 
  # facet_wrap(~Software, labeller = as_labeller(swlabels)) + 
  ylab("Elapsed time (s)") + 
  xlab("Number of genomes")
ggsave("2.4_Analysis_Running-time/Running-time.png")

f2b <- rt_tsvs %>%
  lapply(read.csv, sep="\t", header = TRUE) %>%
  do.call(rbind, .) %>%
  subset(Num_genomes < 150) %>%
  within(Software <- factor(Software, levels = names(software_col))) %>%
  ggplot(aes(x=Num_genomes, y=elapsed_time)) + 
  geom_point(aes(fill=Software), shape = 21, alpha = .3) + 
  stat_summary(fun = mean, geom = "line", mapping = aes(color=Software, group=Software)) +
  scale_color_manual(values = software_col) +
  scale_fill_manual(values = software_col) +
  scale_x_continuous(breaks = NGen) +
  gghighlight::gghighlight() +
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = .5, vjust = .5), 
        strip.background = element_rect(fill = NA), 
        text = element_text(family = "serif")) + 
  facet_wrap(~Software, labeller = as_labeller(swlabels)) + 
  ylab("Elapsed time (s) - Log10 scale") + 
  xlab("Number of genomes")
ggsave("2.4_Analysis_Running-time/Running-time_Log10.png")

library(patchwork)
fig2 <- (f2a + f2b & theme (legend.position = "bottom")) +
  patchwork::plot_annotation(tag_levels = "A") + 
  patchwork::plot_layout(widths = c(1,2), guides = "collect"); fig2
ggsave("2.4_Analysis_Running-time/Figure3.pdf", plot = fig2, width = 25, height = 20, units = "cm")

