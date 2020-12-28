FROM r-base:4.0.3

LABEL maintainer="markus@compbio.cc"
LABEL version="1.0"
LABEL description="Differential Expression pipeline (robust)"

RUN apt update
RUN apt install -y python3 python3-pip git
RUN apt install -y libgit2-dev zlib1g-dev libcurl4-gnutls-dev libxml2-dev libssl-dev libpng-dev libjpeg-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libgl-dev libgsl-dev libcurl4-gnutls-dev libxml2-dev libssl-dev libpng-dev libjpeg-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libgl-dev libgsl-dev 

RUN pip3 install openpyxl jinja2 matplotlib venn pandas numpy seaborn scikit-learn umap scipy statsmodels upsetplot


RUN R -e 'install.packages(c("BiocManager", "devtools", "argparse", "dbplyr"))'
RUN R -e 'BiocManager::install("clusterProfiler")'
RUN R -e 'BiocManager::install("annotables")'
RUN R -e 'BiocManager::install("RDAVIDWebService")'
RUN R -e 'BiocManager::install("qvalue")'
RUN R -e 'BiocManager::install("org.Mm.eg.db")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'BiocManager::install("org.Sc.sgd.db")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("EnrichmentBrowser")'
RUN R -e 'BiocManager::install("statmod")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("data.table")'
RUN R -e 'devtools::install_github("zimmerlab/MS-EmpiRe")'


WORKDIR /git
RUN /bin/bash -c "chmod -R 775 /git"
RUN git clone https://github.com/mjoppich/porestat /git/poreSTAT

ENTRYPOINT ["python3", "/git/poreSTAT/porestat/DEtools/DifferentialAnalysis.py"]
CMD ["python3", "/git/poreSTAT/porestat/DEtools/DifferentialAnalysis.py", "--help"]