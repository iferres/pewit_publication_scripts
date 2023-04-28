library(magrittr)

orgs <- readLines("figshare_data/pangenome/gene_presence_absence.Rtab", n=1) %>% strsplit("\t") %>% extract2(1)
orgs <- orgs[-1]

gffs <- list.files(path="gffs/", pattern="[.]gff$", full.names=T)
gffs <- gffs[sub("[.]gff$","", sub("#", "_", basename(gffs))) %in% orgs] 

hmm <- "Pfam/Pfam-A.hmm"
dat <- "Pfam/Pfam-A.hmm.dat"

library(pewit)

p <- pangenome(gffs, hmm, dat, n_threads = 20)
p$save_pangenomeRDS("pangenome_Serratia2.RDS")
# p <- load_pangenomeRDS("pangenome_Serratia2.RDS")

p$core_level <- 99

library(parallel)
# 1. Align translation of core genes
ali <- p$core_seqs_4_phylo() %>%           # Core genome sequences
	  mclapply(DECIPHER::AlignSeqs, 
		   verbose=FALSE, 
		   mc.cores=20)  # Align


# 3. Concatenate neutral core clusters
concat_ali <- ali %>%            # Select neutral clusters
	do.call(Biostrings::xscat, .) %>%       # Concatenate alignments
	setNames(p$organisms$org) %>%           # Set sequence names
	as('matrix') %>%                        # Transform to matrix
	tolower()                               # Translate to lower case

# Keep only polymorphic sites
concat_ali <- concat_ali[, !apply(concat_ali, 2, function(x) all(x == x[1]))]

saveRDS(concat_ali, "concat_ali.RDS")
# concat_ali <- readRDS("concat_ali.RDS")
writeXStringSet(DNAStringSet(unlist(apply(concat_ali, 1, paste0, collapse="", simplify = F))), filepath = "concat_ali.fasta")

# 4. Compute structure
library(rhierbaps)
rhb <- hierBAPS(snp.matrix = concat_ali, # Input matrix alignment
                n.pops = 100,            # Max number of subpopulations
		max.depth = 4,           # Max depth for hierarchical clustering
		n.extra.rounds = 10,     # Extra rounds to ensure convergence 		
		n.cores=15)

saveRDS(rhb, "rhb.RDS")
rhb <- readRDS("rhb.RDS")

res <- rhb$partition.df
lin <- data.frame(org = as.character(res[, 1]), 
		                    lineage = as.factor(res[, 2]))
p$add_metadata(map = 'org', data = lin)


# 5. Write input data for twilight pipeline
library(magrittr)
p$pan_matrix %>% 
	t() %>% 
	as.data.frame() %>% 
	tibble::rownames_to_column(var="Gene") %>% 
	write.table(file= "gene_presence_absence.Rtab", sep="\t", row.names=FALSE, quote=FALSE)

rhb$partition.df[, 1:2] %>% 
	write.table(file="grouping_rHierbaps_level1.tab", 
		    row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# 6. Download twighlight script and make it executable
download.file(url="https://raw.githubusercontent.com/ghoresh11/twilight/master/classify_genes.R", 
	      destfile="classify_genes.R", 
	      method="wget") 
system("chmod +x classify_genes.R")

# 7. Run twighlight pipeline
system("Rscript ./classify_genes.R -p gene_presence_absence.Rtab -g grouping_rHierbaps_level1.tab")

# Read twilight results and add info to pangenome 
twi <- read.csv("results/out/classification.tab", sep = "\t")
colnames(twi)[1] <- "cluster"
p$add_metadata("cluster", twi)


# 8. Compute ML tree with IQTree, uses the same model as in the paper:
# https://doi.org/10.1038/s41467-022-32929-2
system("singularity exec -B /mnt/cive/nacho pewit_extra.sif iqtree2 -s concat_ali.fasta.varsites.phy -T 15 -m TIM2e+ASC+R4")

# Read tree 
library(ape)
library(ggtree)
library(phangorn)
nj <- read.tree("concat_ali.fasta.varsites.phy.bionj") 
tree <- ggtree(ladderize(midpoint(nj)))

# Add RHierBAPS info to tree
tree <- tree %<+% 
  rhb$partition.df

# Get taxa name order to use it later
txnm <- get_taxa_name(tree)

# Extract panmatrix and convert it to a binary presence/absence matrix
pm <- p$pan_matrix
pm[which(pm >= 1, arr.ind = TRUE)] <- 1L

# The following is to order the columns (clusters) to make the heatmap "pretty"
csm <- colSums(pm)
spl <- split(names(csm), csm)
tpm <- t(pm)
norder <- lapply(spl, function(x) {
  if (length(x)<2){
    x
  }else{
    d <- vegan::vegdist(tpm[x, , drop=F], method = "jaccard", binary = T, na.rm = T)
    hc <- hclust(d, "single")
    x[hc$order]
  }
})
norder <- norder[order(as.integer(names(norder)), decreasing = T)]
forder <- unlist(norder)

# Transform panmatrix data to long format (required to plot with ggplot2)
library(reshape2)
ml <- melt(pm)
colnames(ml)[1:2] <- c("org", "cluster")

# Add twilight classification to the long-formated panmatrix
ml$general_class <- twi$general_class[match(ml$cluster, twi$cluster)]
ml$specific_class <- twi$specific_class[match(ml$cluster, twi$cluster)]

# Factor some data types to manipulate order in figures
ml$value <- as.factor(ifelse(ml$value>=1, "Present", "Absent"))
ml$org <- factor(ml$org, levels = rev(txnm))
ml$cluster <- factor(ml$cluster, levels = forder)
ml$specific_class <- factor(ml$specific_class, 
                            levels = c("Collection core",
                                       "Multi-lineage core",
                                       "Lineage specific core",
                                       "Multi-lineage intermediate",
                                       "Lineage specific intermediate",
                                       "Multi-lineage rare",
                                       "Lineage specific rare",
                                       "Core and intermediate",
                                       "Core and rare",
                                       "Core, intermediate and rare",
                                       "Intermediate and rare", 
                                       "Absent in large lineages"))

# Faceted panmatrix with twilight categories 

library(ggplot2)
ggpa <- ml %>%
  # dplyr::filter(specific_class %in% c("Lineage specific core",
  #                                     "Multi-lineage intermediate")) %>%
  ggplot(aes(x=cluster, y = org, fill = value)) + 
  geom_raster() + 
  scale_fill_grey(start = .9, end = .2)  +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(~specific_class, 
             space = "free", 
             scales = "free_x",
             labeller = labeller(grp = label_wrap_gen(width = 10))) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

# Colour bar showing each twilight category
ggbr <- ml %>%
  # dplyr::filter(specific_class %in% c("Lineage specific core",
  #                                     "Multi-lineage intermediate")) %>%
  ggplot(aes(x = cluster, fill = specific_class)) + 
  geom_bar(stat = "count", width = 1) + 
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(~specific_class, 
             space = "free", 
             scales = "free_x") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        legend.title = element_blank())


# Identify clade nodes
clnd <- sapply(split(rhb$partition.df$Isolate, rhb$partition.df$`level 1`), 
                function(x) MRCA(tree, x))

# Create colors for clades
clade_cols <- colorRampPalette(c(ggsci::pal_npg()(10), ggsci::pal_jco()(10)))(length(clnd))
names(clade_cols) <- names(clnd)
clade_cols <- c(clade_cols, "0"="black")

# Coloured tree
tree2 <- groupClade(tree, clnd, "Lineage") + 
  aes(color=Lineage) + 
  scale_color_manual(values = clade_cols) +
  theme(legend.position = "none")

# Arrange plots
library(patchwork)
gg <- (tree2 +  ggpa + plot_spacer() + ggbr) +
  plot_layout(ncol = 2, 
              widths = c(1, 3), 
              heights = c(20,1),)

# Save
ggsave("raster.png", gg, width = 12)

