sessionInfo()
"
R version 4.0.5 (2021-03-31)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Fedora 34 (Workstation Edition)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libflexiblas.so.3.1

locale:
 [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8    
 [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=es_ES.UTF-8   
 [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.5    CoprManager_0.3.9
"


################
## Containers ##
################

container_dir <- 'Containers'
if (!dir.exists(container_dir)) dir.create(container_dir)

## Repositories ##

container_sources <- list(
  
  gff3toembl = c(transport = "docker",
                 registry = "quay.io",
                 collection = "biocontainers",
                 name = "gff3toembl",
                 tag = "1.1.4--pyh864c0ab_2"),
  
  roary = c(transport = "docker",
            registry = "quay.io",
            collection = "biocontainers",
            name = "roary",
            tag = "3.13.0--pl526h516909a_0"),

  panx = c(transport = "docker",
           registry = "quay.io",
           collection = "biocontainers",
           name = "panx",
           tag = "1.6.0--py27_0"),
  
  panaroo = c(transport = "docker",
              registry = "quay.io",
              collection = "biocontainers",
              name = "panaroo",
              tag = "1.2.10--pyhdfd78af_0"),
  
  pirate = c(transport = "docker",
              registry = "quay.io",
              collection = "biocontainers",
              name = "panaroo",
              tag = "1.2.10--pyhdfd78af_0")
)

# Download 
container_paths <- lapply(container_sources, function(x){
  location <- normalizePath(container_dir)
  subdir <- paste0(location, "/", x["name"])
  if(!dir.exists(subdir)) dir.create(subdir)
  path <- paste0(subdir, "/", x["name"], "_", x["tag"], ".sif")
  if (!file.exists(path)){
    url <- paste0(x["transport"], "://",
                  x["registry"], "/",
                  x["collection"], "/",
                  x["name"], ":",
                  x["tag"])
    run <- paste("singularity pull", path, url)
    system(run) 
  }
  return(path)
})


# Return
container_paths
