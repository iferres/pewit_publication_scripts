

# Compute contingency table function
`%contingency%` <- function(x, y){
  
  #Order
  dn <- dimnames(x)[[1]]
  x <- x[dn, dn]
  y <- y[dn, dn]
  
  #Clean
  diag(x) <- 0L
  diag(y) <- 0L
  
  #Compute contigency
  # tbl <- table(x%checkv%y)
  tbl <- table(factor(paste(x, y, sep = ''), levels = c('00', '01', '10', '11')))
  names(tbl)[which(names(tbl)==c('00', '01', '10', '11'))] <- c('TN','FP','FN','TP')
  #Substract diagonal of "TRUE NEGATIVES"
  dm <- dim(x)
  tbl['TN'] <- tbl['TN'] - dm[2]
  #Values are repeated, so..
  res <- tbl/2
  
  return(res)
}


bench_roary <- function(din, 
                        dout = 'roaryi70_vs_sim', 
                        name = "roary_i70",
                        container = value(containers)$roary,
                        container_exec_args = "-B $PWD",
                        cmd = "roary",
                        cmd_args = "-p 1 -i 70",
                        work_dir = "1.3_Run-clustering-benchmark"){
  
  meta <- strsplit(basename(din), "_")[[1]]
  Ne <- as.numeric(meta[2])
  Replicate <- as.integer(sub("rep", "", meta[3]))
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  f <- paste0(dres, "/", basename(din))  
  #if (dir.exists(f)) unlink(f, recursive = TRUE)
  
  run <- paste(
    "singularity exec",
    container_exec_args,
    container,
    cmd, 
    cmd_args,
    "-f", f, 
    paste(gffs, collapse = ' '), '2> /dev/null'
  )
  
  # Check if already run, avoid re-running if so
  sha1 <- digest::sha1(run)
  sha1fi <- paste0(dres, "/.", basename(din), "_", sha1, ".sha1run")
  if (! file.exists(sha1fi) ) {
    system(run)
    file.create(sha1fi)
  }  
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  
  # Load run
  csv <- paste(f, '/gene_presence_absence.csv', sep = '')
  pg <- pagoo::roary_2_pagoo(csv)
  grp <- lapply(pg$genes, function(x) as.character(x$gene))
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  rm(mat_expect)
  rm(mat_observed)
  gc()
  
  # return
  output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
  return(output)
}


bench_panx <- function(din, 
                       dout = 'panx_vs_sim', 
                       name = "panx",
                       container = value(containers)$panx,
                       container_exec_args = "-B $PWD",
                       cmd = "panX.py",
                       cmd_args = "-sl Simulated_bacteria -cg 0.7 -st 1 2 3 4 5 6 -t 1 -ct",
                       work_dir = "1.3_Run-clustering-benchmark"){
  
  meta <- strsplit(basename(din), "_")[[1]]
  Ne <- as.numeric(meta[2])
  Replicate <- as.integer(sub("rep", "", meta[3]))
  
  gbks <- list.files(path = din, pattern = '[.]gbk$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  fn <- paste0(dres, "/", basename(din))  
  #if (dir.exists(fn)) unlink(fn, recursive = TRUE)
  
  run <- paste(
    "singularity exec",
    container_exec_args,
    container,
    cmd, 
    "-fn", fn,
    cmd_args,
    "> /dev/null"
  )
  
  # Check if already run, avoid re-running if so
  sha1 <- digest::sha1(run)
  sha1fi <- paste0(dres, "/.", basename(din), "_", sha1, ".sha1run")
  if (! file.exists(sha1fi) ) {
    
    dir.create(fn)
    file.copy(normalizePath(gbks), paste(fn, basename(gbks), sep = '/'))
    tr <- try(system(run, ignore.stdout = TRUE, ignore.stderr = TRUE))
    file.create(sha1fi)
    
    if (class(tr)=="try-error"){
      
      res <- list(TN = NA_real_, FP = NA_real_, FN = NA_real_, TP = NA_real_)
      output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
      return(output)
      
    } else {
      
      #clean: remove all files but allclusters_final.tsv and allclusters.tsv
      alls <- list.files(path = fn, recursive = T, full.names = TRUE)
      torm <- grep('allclusters.*[.]tsv$', alls, invert = TRUE, value = TRUE)
      file.remove(torm)
      
      tsvs <- dir(path = fn, recursive = TRUE)
      tsv2 <- paste(fn, basename(tsvs), sep = '/')
      file.rename(paste(fn, tsvs, sep = '/'), tsv2)
      unlink(list.dirs(fn)[-1], recursive = T)
      gpl <- grepl('allclusters_final.tsv', tsv2, fixed = TRUE)
      if (any(gpl)){
        allclusters_final_tsv <- tsv2[which(gpl)]
      }else{
        # stop('No allclusters_final.tsv file.')
        allclusters_final_tsv <- tsv2[1]
      }
      
    }
  
  } else {
    
    allclusters_final_tsv <- paste0(fn, "/allclusters_final.tsv")
    if (!file.exists(allclusters_final_tsv)) {
      res <- list(TN = NA_real_, FP = NA_real_, FN = NA_real_, TP = NA_real_)
      output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
      return(output)
    }
    
  }  
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  panx_2_pagoo <- function(allclusters_final_tsv, group_prefix = 'group', sep = '__'){
    
    x <- readLines(allclusters_final_tsv)
    spl <- strsplit(x, '\t')
    spl <- setNames(spl, paste0(group_prefix, seq_len(length(x))))
    sb <- lapply(spl, function(i) sub('[|]', sep, i))
    df <- reshape2::melt(sb)
    df$org <- sub(paste0(sep, '.+$'), '', df$value)
    colnames(df) <- c('gene', 'cluster', 'org')
    pg <- pagoo::PgR6M$new(data = df, sep = sep)
    return(pg)
  }
  
  pg <- panx_2_pagoo(allclusters_final_tsv = allclusters_final_tsv, sep = '_')
  
  grp <- lapply(pg$genes, function(x) as.character(x$gene))
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  
  # return
  output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
  return(output)
  
}




bench_micropan_blast <- function(din, 
                                 dout = 'micropanBlastComplete_vs_sim', 
                                 name = "micropanBlast_Complete",
                                 linkage = "complete",
                                 job = 1,
                                 work_dir = "1.3_Run-clustering-benchmark"){
  
  meta <- strsplit(basename(din), "_")[[1]]
  Ne <- as.numeric(meta[2])
  Replicate <- as.integer(sub("rep", "", meta[3]))
  
  faas <- list.files(path = din, pattern = '[.]faa$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  dd <- paste0(dres, '/', basename(din))
  dir.create(dd)
  
  
  gid <- paste0("%0",nchar(length(faas)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(faas)))
  ref <- cbind(basename(faas), gid)
  write.csv(ref, file = paste0(dd, '/ref_gid.csv'),quote = F)
  
  prep_faas <- sapply(faas, function(x){
    ff <- basename(x)
    gd <- ref[which(ref[, 1] == ff), 2]
    out <- paste(dd, ff, sep = '/')
    micropan::panPrep(in.file = x, genome_id = gd, out.file = out)
  })
  
  nams <- lapply(prep_faas, function(x){
    fa <- read.fasta(x)
    an <- sapply(fa, attr, 'Annot')
    spl <- strsplit(sub('^>', '', an), ' ')
    mm <- do.call(rbind, spl)
    rownames(mm) <- NULL
    return(mm)
  })
  
  nams <- do.call(rbind, nams)
  
  blout_dir <- paste(dd, 'blast_out', sep = '/')
  # if (dir.exists(blout_dir)) unlink(blout_dir, recursive = TRUE)
  dir.create(blout_dir)
  
  runMicropanBlast <- function(prep_faas, oblf, job = 1, linkage = "complete"){
    micropan::blastpAllAll(prep_faas, 
                          oblf, 
                          e.value = 1e-10,
                          threads = 1L,
                          job = job,
                          verbose = FALSE)
    bl <- list.files(path = oblf, full.names = TRUE)
    df <- micropan::bDist(bl, verbose = FALSE)
    # unlink(oblf, recursive = TRUE)
    bcl <- micropan::bClust(df, threshold = 0.75, linkage = linkage, verbose = FALSE)
    return(bcl)
  }
  
  clust <- runMicropanBlast(prep_faas = prep_faas, 
                            oblf = blout_dir, 
                            job = job, 
                            linkage = linkage)
  
  names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
  grp <- split(names(clust), clust)
  
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  res <- mat_expect %contingency% mat_observed
  # return
  output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
  return(output)
  
}




bench_micropan_hmmer <- function(din, dout = 'micropanHmmer_vs_sim', hmm_pfam){
  
  require(micropan)
  require(seqinr)
  
  faas <- list.files(path = din, pattern = '[.]faa$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  dd <- paste0(dout, '/', din)
  if (dir.exists(dd)) unlink(dd, recursive = TRUE)
  dir.create(dd)
  
  
  gid <- paste0("%0",nchar(length(faas)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(faas)))
  ref <- cbind(basename(faas), gid)
  write.csv(ref, file = paste0(dd, '/ref_gid.csv'),quote = F)
  
  prep_faas <- sapply(faas, function(x){
    ff <- basename(x)
    gd <- ref[which(ref[, 1] == ff), 2]
    out <- paste(dd, ff, sep = '/')
    panPrep(in.file = x, GID.tag = gd, out.file = out)
  })
  
  nams <- lapply(prep_faas, function(x){
    fa <- read.fasta(x)
    an <- sapply(fa, attr, 'Annot')
    spl <- strsplit(sub('^>', '', an), ' ')
    mm <- do.call(rbind, spl)
    rownames(mm) <- NULL
    mm
  })
  
  nams <- do.call(rbind, nams)
  
  hmmout_dir <- paste(dd, 'hmmer_out', sep = '/')
  if (dir.exists(hmmout_dir)) unlink(hmmout_dir, recursive = TRUE)
  dir.create(hmmout_dir)
  
  runMicropanHmmer <- function(prep_faas, oblf, db){
    hmmerScan(in.files = prep_faas, 
              db = db, 
              out.folder = oblf, 
              threads = 4L,
              verbose = FALSE)
    pfam.files <- list.files(path = oblf, full.names = TRUE)
    pfam.table <- lapply(pfam.files, function(i){
      tab <- readHmmer(i)
      tab <- hmmerCleanOverlap(tab)
      tab
    })
    pfam.table <- do.call(rbind, pfam.table)
    
    cluster.pfam <- dClust(pfam.table)
    return(cluster.pfam)
  }
  
  clust <- runMicropanHmmer(prep_faas = prep_faas, oblf = hmmout_dir, db = hmm_pfam)
  names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
  grp <- split(names(clust), clust)
  
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  res <- mat_expect %contingency% mat_observed
  return(res)
  
}



bench_panaroo <- function(din, 
                          dout = "panarooStrict_vs_sim", 
                          name = "panarooStrict",
                          container = value(containers)$panaroo,
                          container_exec_args = "-B $PWD",
                          cmd = "panaroo",
                          cmd_args = "--clean-mode strict --threads 1",
                          work_dir = "1.3_Run-clustering-benchmark",
                          dry_run = FALSE){
  
  meta <- strsplit(basename(din), "_")[[1]]
  Ne <- as.numeric(meta[2])
  Replicate <- as.integer(sub("rep", "", meta[3]))
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  f <- paste0(dres, "/", basename(din))  
  #if (dir.exists(f)) unlink(f, recursive = TRUE)
  
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
  
  if (dry_run) {
    # print(run)
    return(run)
  }
  
  # Check if already run, avoid re-running if so
  sha1 <- digest::sha1(run)
  sha1fi <- paste0(dres, "/.", basename(din), "_", sha1, ".sha1run")
  if (! file.exists(sha1fi) ) {
    system(run)
    file.create(sha1fi)
  }  
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  
  # Load run
  csv <- paste(f, '/gene_presence_absence.csv', sep = '')
  pg <- pagoo::panaroo_2_pagoo(csv)
  grp <- lapply(pg$genes, function(x) as.character(x$gene))
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  rm(mat_expect)
  rm(mat_observed)
  gc()
  
  # return
  output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
  return(output)
}





bench_pewit <- function(din, 
                        dout = "pewit_vs_sim", 
                        name = "pewit",
                        hmm,
                        dat,
                        work_dir = "1.3_Run-clustering-benchmark"){
  
  require(pewit)
  
  meta <- strsplit(basename(din), "_")[[1]]
  Ne <- as.numeric(meta[2])
  Replicate <- as.integer(sub("rep", "", meta[3]))
  
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(work_dir)) dir.create(work_dir)
  
  dres <- paste0(work_dir, "/", dout)
  if (!dir.exists(dres)) dir.create(dres)
  
  dd <- paste0(dres, '/', basename(din))
  dir.create(dd)
  
  fout <- paste0(dd, "/pewit_out_", basename(din), ".RDS")
  
  if (file.exists(fout)){
    
    pg <- pagoo::load_pangenomeRDS(fout)
      
  } else {
    
    if (missing(hmm) | missing(dat)){
      pg <- pangenome(gffs = gffs, 
                      n_threads = 1, 
                      verbose = F)
    } else {
      pg <- pangenome(gffs = gffs, 
                      hmm_pfam = hmm, 
                      dat_pfam = dat, 
                      n_threads = 1, 
                      verbose = F)
    }
    pg$save_pangenomeRDS(file = fout)
  }
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  grp <- lapply(pg$genes, function(x) as.character(x$gene))
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  
  # return
  output <- as.data.frame(c(as.list(res), "Software" = name, "Ne"=Ne, "Replicate" = Replicate))
  return(output)
  
}




bench_pirate <- function(din, 
                         dout = 'pirate_vs_sim', 
                         cmd = 'singularity run singularity_containers/pirate.sif',
                         cmd2 = 'singularity exec singularity_containers/pirate.sif PIRATE_to_roary.pl'){
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  f <- paste0(dout, '/', din)
  if (dir.exists(f)) unlink(f, recursive = TRUE)
  run <- paste(cmd, '-k "--diamond" --para-off -t 4 -o', f,'-i', din)
  system(run)
  run2 <- paste0(cmd2, " -i ", f, "/PIRATE.*.tsv -o ", f, "/gene_presence_absence.csv")
  system(run2)
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  
  # Load run
  csv <- paste(f, '/gene_presence_absence.csv', sep = '')
  pg <- pagoo::roary_2_pagoo(csv, paralog_sep = ";")
  df <- unlist(pg$gene, use.names = FALSE)
  grp <- split(as.character(df$gene), df$group)
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  return(res)
}


