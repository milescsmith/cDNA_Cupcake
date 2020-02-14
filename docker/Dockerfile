FROM python:3.8-slim-buster

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y \
        gcc \
        git && \
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/* && \
    rm -rf /var/lib/apt/lists/*

RUN pip --no-cache-dir install cython numpy && \
    pip --no-cache-dir install biopython scikit-learn pysam bx-python && \
    pip --no-cache-dir install git+https://github.com/milescsmith/cDNA_Cupcake
