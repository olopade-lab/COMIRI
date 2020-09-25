### Modified from jemimalwh/trust4 ###

FROM ubuntu:16.04

ARG TRUST4_VERSION=v1.0.0

ENV PACKAGES apt-utils git

ENV TRUST4PATH /home/TRUST4

RUN set -ex

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    cd /home && \
    export GIT_SSL_NO_VERIFY=1 && \
    git clone -b ${TRUST4_VERSION} --single-branch https://github.com/liulab-dfci/TRUST4.git && \
    mkdir /work && cd /work && \
    apt-get --purge autoremove -y ${PACKAGES}

WORKDIR /work

ENV PATH /home/TRUST4:${PATH}