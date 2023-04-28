

ffn2gff <- function(ffn){
  require(seqinr)
  rf <- read.fasta(ffn, 
                   seqtype = 'DNA', 
                   as.string = TRUE, 
                   forceDNAtolower = TRUE)
  
  li1 <- '##gff-version 3\n'
  
  contig <- 1:length(rf)
  len_sq <- sapply(rf, nchar)
  li2 <- paste0('##sequence-region contig.', contig, ' 1 ', len_sq, '\n')
  
  nn <- names(rf)
  li3 <- paste0('contig.', contig, 
                '\tProdigal:2.6\tCDS\t1\t',
                len_sq,
                '\t.\t+\t0\tID=',nn, 
                ';locus_tag=',nn, 
                ';gene=',nn,
                ';product=',nn,'\n') 
  li4 <- '##FASTA\n'
  li5 <- paste0('>contig.',contig,'\n',unlist(rf),'\n')
  
  o <- sub('[.]f\\w+$', '.gff', ffn)
  cat(c(li1, li2, li3, li4, li5), file = o, sep = '')
  
  return(o)
}


gff3_2_embl <- function(gff3, 
                        container = value(containers)$gff3toembl,
                        cmd = 'gff3_to_embl',
                        authors = 'John', 
                        title = 'Some Title',
                        publication = 'Some Journal', 
                        genome_type = 'circular', 
                        classification = 'PROK',
                        organism = 'Organism', 
                        taxonid = 1234, 
                        accession = 'project',
                        description = 'description'){
  
  arg <- paste0("--authors ", "\'",authors,"\'", ' ',
                "--title ", "\'", title, "\'", ' ',
                "--publication ", "\'", publication, "\'", ' ',
                "--genome_type ","\'" , genome_type, "\'", ' ',
                "--classification ","\'" , classification, "\'", ' ',
                "--output_filename /dev/stdout", ' ',
                organism,  ' ',
                taxonid,  ' ',
                "\'", accession, "\'",  ' ',
                "\'", description, "\'", ' ',
                gff3)
  
  run <- paste("singularity exec -B $PWD", container, cmd, arg)
  
  stout <- system(run, intern = T)
  
  deli <- grep('//', stout)
  
  fctr <- rep(seq_along(deli), c(deli[1], diff(deli)))
  
  chunks <- split(stout, fctr)
  
  # fields to add 
  DT1 <- 'DT   27-MAR-2019 (Rel. 01, Created)'
  DT2 <- 'DT   27-MAR-2019 (Rel. 02, Last updated, Version 1)'
  genomeid <- sub('[.]gff$', '', basename(gff3))
  OS <- paste0('OS   Fakeum bacterii (', genomeid, ')')
  
  chunks <- lapply(chunks, function(x){
    
    ln <- length(x)
    out <- vector('character', ln + 4 + 3)
    gpac <- rev(grep('^AC', x))[1]
    out[1:(gpac+1)] <- x[1:(gpac+1)]
    out[(gpac+2):(gpac+4)] <- c(DT1, DT2, 'XX')
    gpde <- rev(grep('^DE', x))[1]
    rg <- (gpac+2):(gpde+1)
    lim <- (gpac+5 + length(rg)-1)
    out[(gpac+5):lim] <- x[rg]
    ini <- lim+1
    lim <- ini + 1
    KW <- 'KW   function'
    out[ini:lim] <- c(KW, 'XX')
    ini <- lim+1
    lim <- ini+1
    out[ini:lim] <- c(OS, 'XX')
    ini <- lim+1
    lim <- ln + 4 + 3
    out[ini:lim] <- x[(gpde+2):length(x)]
    out
    
  })
  
  
  
  embl <- sub('[.]gff$', '.embl', gff3)
  if (file.exists(embl)){
    file.remove(embl)
  }
  cat(paste0('Writing embl at ', embl, '...'))
  lapply(chunks, cat, file = embl, sep = '\n', append=TRUE)
  if (file.exists(embl)) cat(' DONE!\n') else cat(' FAIL :(\n')
  embl
}



embl_2_gbk <- function(embl){
  
  rlembl <- readLines(embl)
  sp <- grep('^//$', rlembl)
  lsp <- length(sp)
  ini <- c(1L, (sp+1)[-lsp])
  end <- sp
  mm <- cbind(ini, end)
  chunks <- apply(mm, 1, function(x){rlembl[x[1]:x[2]]})
  
  # if (missing(file)) file <- stdout()
  #Line 1
  
  process_chunk <- function(x){
    contig <- sub('[^contig\\d+]+', '', grep('^AC.+contig\\d+', x, value = TRUE))
    large <- regmatches(x[1], regexpr('\\d+ BP.', x[1]))
    large <- sub('BP[.]', 'bp', large)
    L1 <- paste('LOCUS       ', 
                contig, '             ',
                large, '    ',
                'DNA     ', 
                'circular ', 
                '01-APR-2019', 
                sep='')
    
    #Line 2
    L2 <- paste('DEFINITION  Genus species strain strain.')
    
    L3 <- paste('ACCESSION   ')
    
    L4 <- paste('VERSION     ')
    
    #Line 5
    KW <- sub('^KW[ ]+','',grep('^KW[ ]+', x, value = TRUE))
    L5 <- paste('KEYWORDS    ', KW, sep = '')
    
    #Line 6
    organism <- sub('^OS[ ]+','',grep('^OS[ ]+', x, value = TRUE))
    L6 <- paste('SOURCE      ', organism, sep = '')
    
    L7 <- paste('  ORGANISM  ', organism, sep = '')
    
    L8 <- paste('COMMENT     ', 'This is a fake GenBank file.', sep='')
    
    features <- gsub('FH|[ ]|Key','',grep('^FH[ ]+Key[ ]+', x, value = TRUE))
    L9 <- paste('FEATURES             ', features, sep = '')
    
    L10 <- sub('^FT','  ',grep('^FT', x, value = T))
    
    
    #L11: translation 
    
    sq <- grep('^SQ', x)
    L12 <- 'ORIGIN'
    
    dna <- x[(sq+1):(length(x)-1)]
    string <- paste0(gsub('[ ]|\\d', '', dna), collapse = '')
    trans <- paste0(seqinr::translate(strsplit(string, '')[[1]]), collapse = '')
    trans <- sub('[*]', '', trans)
    translation <- paste0('/translation="', trans, '"')
    splitInParts <- function(string, size){
      pat <- paste0('(?<=.{',size,'})')
      strsplit(string, pat, perl=TRUE)[[1]]
    }
    translation <- splitInParts(translation, 59)
    space <- sub('[^ ]+', '', L10)[2]
    L11 <- paste0(space, translation, sep = '')
    
    counts <- cumsum(c(1, rep(60, length(dna)-1)))
    mx <- max(nchar(counts))
    dna <- sub('[ ]+\\d+', '', dna)
    dna <- sub('^[ ]+','',dna)
    spp <- sapply((9-nchar(counts)), function(x){
      paste0(rep(' ', x), collapse = '')
    })
    L13 <- paste(paste(spp, counts, ' ', sep=''), dna, sep='')
    L14 <- x[length(x)]
    
    tofile <- c(L1, L2, L3, L4,
                L5, L6, L7, L8, 
                L9, L10, L11, L12, 
                L13, L14)
    tofile
    
    # out <- sub('[.]embl$', '.gbk', embl)
    # cat(tofile, file = out,  sep = '\n', append = TRUE)
  }
  
  out <- sub('[.]embl$', '.gbk', embl)
  tofile <- lapply(chunks, process_chunk)
  if (file.exists(out)) file.remove(out)
  lapply(tofile, function(x) cat(x, file = out, sep = '\n', append = TRUE))
  out
}
