### Modified from jemimalwh/trust4 ###

FROM ubuntu:16.04

ARG TRUST4_VERSION=v1.0.0

RUN apt-get update && \
    apt-get install -y build-essential \
    libpthread-stubs0-dev git perl zlib1g-dev && \
    apt-get clean && \
    export GIT_SSL_NO_VERIFY=1 && \
    git clone -b ${TRUST4_VERSION} --single-branch https://github.com/liulab-dfci/TRUST4.git && \
    cd /TRUST4 && make

ENV PATH "$PATH:/TRUST4"

