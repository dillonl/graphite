FROM ubuntu:latest
MAINTAINER dillonl https://github.com/dillonl

LABEL vender="dillonl"

RUN apt-get update && apt-get -y install \
		cmake \
		git \
		gcc \
		build-essential \
		python-dev \
		libpng-dev
        
## Pull repo
RUN git clone https://github.com/dillonl/graphite.git
RUN mkdir graphite/bin
RUN cd graphite/bin && cmake ../
RUN cd graphite/bin && make

