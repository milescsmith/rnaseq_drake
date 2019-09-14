# FROM rocker/r-ver:3.6.1
# LABEL maintainer="miles-smith@omrf.org"
# RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
#   && apt-get install -y git-core \
# 	libcurl4-openssl-dev \
# 	libpng-dev \
# 	libssl-dev \
# 	libxml2-dev \
# 	make \
# 	pandoc \
# 	pandoc-citeproc \
# 	python-minimal && \
#   apt-get clean && \
#   rm -rf /tmp/downloaded_packages/* && \
#   rm -rf /var/lib/apt/lists/*
# RUN install2.r \
# 	tidyverse rlang magrittr pheatmap kableExtra IRdisplay janitor remotes \
# 	RColorBrewer paletteer cowplot irlba ggforce uwot formattable drake rsvd
# RUN install2.r \
# 	-r https://bioconductor.org/packages/3.9/bioc \
# 	-r https://bioconductor.org/packages/3.9/data/annotation \
# 	-r https://bioconductor.org/packages/3.9/data/experiment \
# 	-r https://bioconductor.org/packages/3.9/workflows \
# 	BiocParallel \
# 	clusterProfiler \
# 	DESeq2 \
# 	S4Vectors \
# 	SummarizedExperiment \
# 	tximport
# RUN installGithub.r \
# 	atlanhq/flyio@cea5f645ce66bc5fb60ac7f7f2840fc8593da02e \
# 	milescsmith/moduleScoreR
# WORKDIR /payload/
# CMD [R]

FROM registry.gitlab.com/milothepsychic/r-alpine:latest
LABEL maintainer="miles-smith@omrf.org"

RUN apk --no-cache --update-cache add \
	libxml2-dev \
	git \
	freetype-dev 

RUN Rscript -e "install.packages(c('tidyverse', 'rlang', 'magrittr', 'pheatmap', 'kableExtra', 'IRdisplay', 'janitor', 'remotes', 'RColorBrewer', 'paletteer', 'cowplot', 'irlba', 'ggforce', 'uwot','formattable', 'drake', 'rsvd', 'BiocManager'))" && \
  Rscript -e "BiocManager::install(c('BiocParallel', 'clusterProfiler', 'DESeq2', 'S4Vectors', 'SummarizedExperiment', 'tximport'))"

RUN Rscript -e "BiocManager::install(c('atlanhq/flyio@cea5f645ce66bc5fb60ac7f7f2840fc8593da02e', 'milescsmith/moduleScoreR'))"

CMD ["Rscript"]