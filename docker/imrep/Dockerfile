FROM python:2.7.18

RUN export GIT_SSL_NO_VERIFY=1 && \
    git clone https://github.com/mandricigor/imrep.git && \
    cd /imrep && \
    pip install pysam && \
    pip install biopython==1.76 && \
    pip install intervaltree && \
    pip install jellyfish && \
    pip install numpy && \
    pip install networkx && \
    tar -xf suffix_tree-2.1.tar.gz && \
    cd suffix_tree-2.1 && \
    python setup.py install && \
    cd .. && \
    rm -rf suffix_tree-2.1