# Mapping and Sort

FROM ubuntu:latest
MAINTAINER pzweuj

# aliyun apt source (not necessary)
COPY aliyun.txt /etc/apt/sources.list

# apt
RUN apt-get update && \
    apt-get install -y \
    zlib1g-dev \
    libbz2-dev \
    wget \
    build-essential \
    libncursesw5-dev \
    liblzma-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/list/* /tmp/* /var/tmp/*


WORKDIR /software

# samtools 1.12
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar -jxvf samtools-1.12.tar.bz2 && \
    rm samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix $(pwd) && \
    make && \
    cp samtools /usr/local/bin/

# bwa 0.7.17
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -jxvf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && \
    rm -rf /software/bwa-0.7.17 && \
    rm -rf /software/bwa-0.7.17.tar.bz2

# sambamba 0.8.0
RUN wget https://github.com/biod/sambamba/releases/download/v0.8.0/sambamba-0.8.0-linux-amd64-static.gz && \
    gunzip sambamba-0.8.0-linux-amd64-static.gz && \
    mv sambamba-0.8.0-linux-amd64-static /usr/local/bin/ && \
    chmod 777 /usr/local/bin/sambamba-0.8.0-linux-amd64-static && \
    ln -s /usr/local/bin/sambamba-0.8.0-linux-amd64-static /usr/local/bin/sambamba && \
    rm -rf /software

# defualt command
WORKDIR /usr/local/bin/
ENV PATH="/usr/local/bin/:${PATH}"
CMD ["/bin/bash"]

