Bootstrap: docker
From: rocker/r-devel

%post
	apt update -y && apt upgrade -y
	apt install -y build-essential libssl-dev hmmer mcl iqtree
	R -e 'install.packages(c("remotes", "rhierbaps", "pegas"))'
	R -e 'remotes::install_github("iferres/pewit")'
