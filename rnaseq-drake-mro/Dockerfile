FROM milescsmith/mro:3.5.3-verse
LABEL maintainer="miles-smith@omrf.org"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y git-core \
	libcurl4-openssl-dev \
	libglu1-mesa-dev \
	libpng-dev \
	libssl-dev \
	libxml2-dev \
	make \
	pandoc \
	pandoc-citeproc \
	python-minimal && \
  apt-get clean -y && \
  apt-get autoremove -y && \
  rm -rf /tmp/downloaded_packages/* && \
  rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('rlang', 'magrittr', 'pheatmap', 'kableExtra', 'IRdisplay', 'janitor', 'remotes', 'RColorBrewer', 'paletteer', 'cowplot', 'irlba', 'ggforce', 'uwot','formattable', 'drake', 'rsvd'))" && \
	Rscript -e "BiocManager::install(c('BiocParallel', 'clusterProfiler', 'DESeq2', 'S4Vectors', 'SummarizedExperiment', 'tximport', 'atlanhq/flyio', 'milescsmith/moduleScoreR', 'tidyverse'))"

WORKDIR /root

CMD ["Rscript"]