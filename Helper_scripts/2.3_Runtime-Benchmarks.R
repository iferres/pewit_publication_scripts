runtime_pewit <- function(gffs, 
                          hmm_pfam, 
                          dat_pfam){
  
  st <- system.time(pg <- try(pewit::pangenome(gffs = gffs,
                                               hmm_pfam = hmm_pfam,
                                               dat_pfam = dat_pfam, 
                                               n_threads = 1L, 
                                               verbose = FALSE)))
  cls <- class(pg)
  exit_code <- ifelse('PgR6MS' %in% cls, 0L, 1L)
  list(system_time = st, exit = exit_code)
}



runtime_roary <- function(gffs,
                          set_id,
                          container = value(containers)$roary,
                          container_exec_args = "-B $PWD",
                          cmd = "roary",
                          cmd_args = "-p 1 -i 70",
                          dout = "roary_i70",
                          work_dir = "2.3_Runtime-Benchmarks"){
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  f <- paste0(dres, '/', set_id)
  if (dir.exists(f)) unlink(f, recursive = TRUE)
  
  run <- paste(
    "singularity exec",
    container_exec_args,
    container,
    cmd, 
    cmd_args,
    "-f", f, 
    paste(gffs, collapse = " "), "2> /dev/null"
  )
  
  st <- system.time(system(run))
  out_csv <- paste(f, "/gene_presence_absence.csv", sep = "") 
  
  exit_code <- ifelse(file.exists(out_csv), 0L, 1L)
  
  list(system_time = st, exit = exit_code)
}


runtime_panx <- function(gbks,
                         set_id,
                         container = value(containers)$panx,
                         container_exec_args = "-B $PWD",
                         cmd = "panX.py",
                         cmd_args = "-sl Cfetus -cg 0.5 -st 1 2 3 4 5 6 -t 1 -ct",
                         dout = "panx",
                         work_dir = "2.3_Runtime-Benchmarks"
                         ){
  
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  fn <- paste0(dres, '/', set_id)
  if (dir.exists(fn)) unlink(fn, recursive = TRUE)
  dir.create(fn)
  
  file.copy(normalizePath(gbks), paste(fn, basename(gbks), sep = '/'))
  
  run <- paste(
    "singularity exec",
    container_exec_args,
    container,
    cmd, 
    "-fn", fn,
    cmd_args,
    "> /dev/null"
  )
  
  st <- system.time(tr <- try(system(run, intern = F, ignore.stdout = T, ignore.stderr = T)))
  
  allclusters <- list.files(path = fn, recursive = TRUE, pattern = "allclusters.*[.]tsv")
  exit_code <- as.integer(length(allclusters)==0)
  
  list(system_time = st, exit = exit_code)
}


runtime_panaroo <- function(
    gffs,
    set_id,
    container = value(containers)$panaroo,
    container_exec_args = "-B $PWD",
    cmd = "panaroo",
    cmd_args = "--clean-mode strict --threads 1",
    dout = "panarooStrict",
    work_dir = "2.3_Runtime-Benchmarks"
){
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  f <- paste0(dres, '/', set_id)
  if (dir.exists(f)) unlink(f, recursive = TRUE)
  
  run <- paste(
    "singularity exec",
    container_exec_args,
    container,
    cmd, 
    cmd_args,
    "-o", f, 
    "-i", paste(gffs, collapse = ' '), 
    "--quiet"
  )
  
  st <- system.time(system(run))
  out_csv <- paste(f, "/gene_presence_absence.csv", sep = "") 
  
  exit_code <- ifelse(file.exists(out_csv), 0L, 1L)
  
  list(system_time = st, exit = exit_code)
  
}


runtime_Micropan_BlastAllAll <- function(
    faas,
    set_id,
    dout = "micropan_BlastAllAll",
    work_dir = "2.3_Runtime-Benchmarks"
  ){
  
  # require(micropan)
  # require(seqinr)
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  fn <- paste0(dres, '/', set_id)
  if (dir.exists(fn)) unlink(fn, recursive = TRUE)
  dir.create(fn)
  
  gid <- paste0("%0",nchar(length(faas)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(faas)))
  ref <- cbind(basename(faas), gid)
  write.csv(ref, file = paste0(fn, '/ref_gid.csv'),quote = F)
  
  prep_faas <- sapply(faas, function(x){
    ff <- basename(x)
    gd <- ref[which(ref[, 1] == ff), 2]
    out <- paste(fn, ff, sep = '/')
    micropan::panPrep(in.file = x, genome_id = gd, out.file = out)
  })
  
  
  blout_dir <- paste(fn, 'blast_out', sep = '/')
  if (dir.exists(blout_dir)) unlink(blout_dir, recursive = TRUE)
  dir.create(blout_dir)
  
  st <- system.time({
    mp <- micropan::blastpAllAll(prep_faas, 
                blout_dir, 
                e.value = 1e-10,
                threads = 1L,
                verbose = FALSE)
  })
  
  exit_code <- as.integer(mp != TRUE)
  
  return(list(system_time = st, exit = exit_code))
}


runtime_Micropan_Clust <- function(
    set_id,
    linkage = "single",
    din = "micropan_BlastAllAll",
    work_dir = "2.3_Runtime-Benchmarks"
  ){
  
  fn <- paste0(work_dir, "/", din, "/", set_id)
  blout_dir <- paste(fn, "blast_out", sep = '/')
  
  bl <- list.files(path = blout_dir, full.names = TRUE)
  st <- system.time({
    df <- micropan::bDist(blast.files = bl, verbose = FALSE)
    bcl <- try(micropan::bClust(df, threshold = 0.75, linkage = linkage, verbose = FALSE))
  })
  
  blastAllAllFile <- paste0(work_dir, "/", "runtime_", din, ".tsv")
  blastAllAllTable <- read.table(blastAllAllFile, header = T)
  st_blast <- blastAllAllTable$elapsed_time[which(blastAllAllTable$set_id == set_id)]
  
  st <- st + st_blast
  list(system_time = st, exit = ifelse(class(bcl)=="try-error", 1L, 0L))
}






