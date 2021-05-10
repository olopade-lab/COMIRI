FROM continuumio/miniconda2:4.7.12

RUN apt-get update && \
    apt-get -y install gcc && \
    conda install -c bioconda bowtie2 && \
    export GIT_SSL_NO_VERIFY=1 && \
    git clone https://github.com/uubram/RTCR.git && \
    cd /RTCR && \
    python setup.py install && \
    rtcr Config Aligner.location=/opt/conda/bin/