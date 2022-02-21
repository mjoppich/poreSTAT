FROM ubuntu:20.04

LABEL maintainer="markus@compbio.cc"
LABEL version="1.5"
LABEL description="Differential Expression pipeline (robust)"

RUN apt update && apt install -y dirmngr apt-transport-https ca-certificates gnupg2
RUN DEBIAN_FRONTEND="noninteractive" apt -y install tzdata

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN apt update && apt install -y python3-dev python3 python3-pip curl git cmake r-base

RUN apt install -y python3-llvmlite libgit2-dev zlib1g-dev libcurl4-gnutls-dev libxml2-dev libhts-dev libssl-dev libpng-dev libjpeg-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libgl-dev libgsl-dev libcurl4-gnutls-dev libxml2-dev libssl-dev libpng-dev libjpeg-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libgl-dev libgsl-dev 
RUN pip3 install matplotlib==3.3.2 numpy
RUN pip3 install umap-learn 
RUN pip3 install matplotlib_venn adjustText natsort openpyxl XlsxWriter jinja2 biopython venn pandas numpy seaborn scikit-learn scipy statsmodels upsetplot HTSeq pysam dill pathos openpyxl h5py

RUN curl -fsSL https://deb.nodesource.com/setup_17.x | bash -
RUN apt-get install -y nodejs

WORKDIR /git
RUN /bin/bash -c "chmod -R 775 /git"

RUN git clone https://github.com/mjoppich/mpld3 /git/mpld3
RUN cd /git/mpld3 && npm update
RUN cd /git/mpld3 && python3 setup.py submodule
RUN cd /git/mpld3 && python3 setup.py buildjs
RUN cd /git/mpld3 && python3 setup.py build
RUN cd /git/mpld3 && python3 setup.py install
RUN cp -r /git/mpld3/mplexporter/mplexporter/renderers/ /usr/local/lib/python3.8/dist-packages/mpld3-0.3.1.dev1-py3.8.egg/mpld3/mplexporter/ && cd /git

RUN git clone https://github.com/mjoppich/porestat /git/poreSTAT
RUN mkdir poreSTAT/clib/build && cd poreSTAT/clib/build && cmake .. && make && cd /git

RUN R -e 'install.packages(c("BiocManager", "devtools", "argparse", "dbplyr", "ashr", "svglite"))'
RUN R -e 'devtools::install_github("stephenturner/annotables")'
RUN R -e 'options(download.file.method = "curl"); devtools::install_github("YuLab-SMU/ReactomePA")'
RUN R -e 'BiocManager::install("clusterProfiler")'
#RUN R -e 'BiocManager::install("ReactomePA")'
RUN R -e 'BiocManager::install("RDAVIDWebService")'
RUN R -e 'BiocManager::install("qvalue")'
RUN R -e 'BiocManager::install("apeglm")'
RUN R -e 'BiocManager::install("org.Mm.eg.db")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'BiocManager::install("org.Sc.sgd.db")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("EnrichmentBrowser")'
RUN R -e 'BiocManager::install("statmod")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("data.table")'
RUN R -e 'devtools::install_github("zimmerlab/MS-EmpiRe")'

ENTRYPOINT ["python3", "/git/poreSTAT/porestat/DEtools/DifferentialAnalysis.py"]
CMD ["python3", "/git/poreSTAT/porestat/DEtools/DifferentialAnalysis.py", "--help"]

#docker build .
#docker build --no-cache .

# => => writing image sha256:09064d1399ebfdf78b45fe3f20d75fe213343eceeeb2e2bf445c5306ad6f2c53
# docker image ls
#docker tag 09064d1399ebfdf78b45fe3f20d75fe213343eceeeb2e2bf445c5306ad6f2c53 mjoppich/porestat_de:latest
#docker tag 09064d1399ebfdf78b45fe3f20d75fe213343eceeeb2e2bf445c5306ad6f2c53 mjoppich/porestat_de:v1.4
#docker login
#docker push mjoppich/porestat_de:v1.4
#docker push mjoppich/porestat_de:latest