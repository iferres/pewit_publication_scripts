simulate_pangenome <- function(
    
    # Roary's publication ref file
    # Kallonen, Teemu (2017): Pan genome reference fasta. figshare. Dataset.
    # https://doi.org/10.6084/m9.figshare.4873022.v1
    # https://figshare.com/articles/Pan_genome_reference_fasta/4873022
    ref = "Data/Ecoli_pan_genome_reference.fa",
    norg = 10L,
    ne = 2e+09, 
    C = 100,
    u = 1E-8,
    v = 1E-11,
    mu = 5E-12,
    write_by = "genome",
    replace = TRUE,
    verbose = FALSE,
    repli = 1L, 
    work_dir = "1.1_Run-simulations"
    
){
  
  require(magrittr)
  
  if (!file.exists(ref)) stop('You have to manually downlowad pangenome reference file, name it 
       "Ecoli_pan_genome_reference.fa", and place it on the Data directory.')
  
  if (! "simurg" %in% installed.packages()[, 1]) stop("simurg is not installed. 
                                                      Use remotes::install_github('iferres/simurg')")
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  Nexp <- sub('e[+]', 'E', ne)
  id <- paste0(Nexp, '_rep', repli)
  dout <- paste0(work_dir, '/sim_', id)
  pg <- simurg::simpg(ref = ref,
                      norg = norg,
                      ne = ne,
                      C = C,
                      u = u,
                      v = v,
                      mu = mu,
                      write_by = write_by,
                      dir_out = dout,
                      replace = replace,
                      verbose = verbose)
  rds <- paste(dout, '/pg_', id, '.RDS', sep = '')
  saveRDS(pg, file = rds)
  message(rds)
  return(rds)
  
}