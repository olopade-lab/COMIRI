FROM ubuntu:16.04

RUN apt-get update && \
    apt-get install -y build-essential \
    libtool automake zlib1g-dev libbz2-dev pkg-config git && \
    export GIT_SSL_NO_VERIFY=1 && \
    git clone https://github.com/neufeld/pandaseq.git && \
    cd /pandaseq && \
    ./autogen.sh && ./configure && make && make install && ldconfig