FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install spliceAI
RUN pip install tensorflow && \
    pip install spliceai

WORKDIR /data
CMD ["spliceai", "-h"]
