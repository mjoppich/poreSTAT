FROM r-base:4.0.3

LABEL maintainer="markus@compbio.cc"
LABEL version="1.0"
LABEL description="Differential Expression pipeline (robust)"


RUN apt update && apt install -y python3-dev python3 python3-pip git cmake
RUN apt install -y python3-llvmlite libgit2-dev zlib1g-dev libcurl4-gnutls-dev libxml2-dev libhts-dev libssl-dev libpng-dev libjpeg-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libgl-dev libgsl-dev libcurl4-gnutls-dev libxml2-dev libssl-dev libpng-dev libjpeg-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libgl-dev libgsl-dev 
RUN apt install -y llvm-10-dev && LLVM_CONFIG=/usr/bin/llvm-config-10 pip3 install umap-learn
RUN pip3 install matplotlib==3.3.2 matplotlib_venn openpyxl jinja2 biopython umap-learn umap venn pandas numpy seaborn scikit-learn scipy statsmodels upsetplot HTSeq pysam dill pathos openpyxl h5py


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

RUN git clone https://github.com/mjoppich/mpld3 /git/mpld3
RUN cd /git/mpld3 && python3 setup.py submodule && python3 setup.py build && python3 setup.py install && cp -r mplexporter/mplexporter/renderers/ /usr/local/lib/python3.9/dist-packages/mpld3-0.3.1.dev1-py3.9.egg/mpld3/mplexporter/ && cd /git

RUN git clone https://github.com/mjoppich/porestat /git/poreSTAT



RUN mkdir poreSTAT/clib/build && cd poreSTAT/clib/build && cmake .. && make && cd /git

ENTRYPOINT ["python3", "/git/poreSTAT/porestat/DEtools/DifferentialAnalysis.py"]
CMD ["python3", "/git/poreSTAT/porestat/DEtools/DifferentialAnalysis.py", "--help"]

#docker build .
#docker tag d426563eec45 mjoppich/porestat_de:latest
#docker tag d426563eec45 mjoppich/porestat_de:v1.3
#docker login
#docker push mjoppich/porestat_de:v1.3 