FROM milescsmith/mro:mro-rnaseq
LABEL maintainer="miles-smith@omrf.org"
RUN \
  export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
    && apt-get install -y \
    git-core \
    libpng-dev \
    pandoc \
    pandoc-citeproc && \
  apt-get clean -y && \
  apt-get autoremove -y && \
  rm -rf /tmp/downloaded_packages/* && \
  rm -rf /var/lib/apt/lists/*

RUN \
  Rscript -e "BiocManager::install(c( \
    'alastairrushworth/inspectdf', \
    'visNetwork', \
    'atlanhq/flyio', \
    'milescsmith/moduleScoreR', \
    'drake', \
    'formattable'))"

WORKDIR /root

CMD ["Rscript"]