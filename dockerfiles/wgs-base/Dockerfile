FROM ubuntu:xenial

# Last updated: September 19, 2018
MAINTAINER Heather Ward <heather@dnastack.com>

# Install dependencies
RUN apt-get update -qq && \
	apt-get install -qq \
	python3-dev \
	python-minimal \
	wget \
	make \
	gcc \
	build-essential \
	libz-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	default-jre \
	default-jdk \
	git


# Install required programs
# bwa v0.7.17
RUN wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz && \
	tar -zxvf v0.7.17.tar.gz && \
	rm v0.7.17.tar.gz
RUN cd bwa-0.7.17 && \
	make

# samtools v1.8
RUN wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
	tar -jxvf samtools-1.8.tar.bz2 && \
	rm samtools-1.8.tar.bz2
RUN cd samtools-1.8 && \
	./configure --prefix=/usr/local/bin && \
	make && \
	make install

# PICARD v2.18.9
RUN wget -P /usr/java/ https://github.com/broadinstitute/picard/releases/download/2.18.9/picard.jar

# GATK v4.0.8.1
RUN git clone --recursive https://github.com/broadinstitute/gatk.git
RUN wget -P /lfs https://github.com/git-lfs/git-lfs/releases/download/v2.5.2/git-lfs-linux-amd64-v2.5.2.tar.gz && \
	cd /lfs && tar -zxvf /lfs/git-lfs-linux-amd64-v2.5.2.tar.gz && \
	rm /lfs/git-lfs-linux-amd64-v2.5.2.tar.gz && \
	/lfs/install.sh
RUN cd /gatk && ./gradlew localJar

# Sambamba v0.6.7
RUN wget -P /usr/bin https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2 && \
	tar -jxvf /usr/bin/sambamba_v0.6.7_linux.tar.bz2 && \
	rm /usr/bin/sambamba_v0.6.7_linux.tar.bz2

# Install tabix/bgzip
RUN apt install -qq tabix

# Add programs to the PATH
ENV PATH $PATH:/bwa-0.7.17/:/usr/local/bin/bin:/gatk:/
ENV PICARD /usr/java/picard.jar
