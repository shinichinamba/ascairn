FROM ubuntu:24.04
LABEL maintainer="Yuichi Shiraishi <friend1ws@gmail.com>" 

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update && apt install -y \
    git \
    wget \
    bzip2 \
    make \
    cmake \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev \
    python3 \
    python3-pip

RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
	tar jxvf samtools-1.21.tar.bz2 && \
	cd samtools-1.21 && \
	./configure && make && make install && \
	cd .. && rm -rf samtools-1.21 samtools-1.21.tar.bz2

RUN apt update && apt install -y jellyfish

ARG TARGETARCH
RUN set -eux; \
    if [ "$TARGETARCH" = "amd64" ]; then \
        wget https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth && \
        chmod +x mosdepth && mv mosdepth /usr/local/bin/ ; \
    elif [ "$TARGETARCH" = "arm64" ]; then \
        wget -qO /usr/local/bin/micromamba \
            https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-aarch64 && \
        chmod +x /usr/local/bin/micromamba && \
        micromamba create -y -p /opt/conda -c conda-forge -c bioconda \
            "mosdepth>=0.3.10" && \
        ln -s /opt/conda/bin/mosdepth /usr/local/bin/mosdepth ; \
    else \
        echo "unsupported TARGETARCH: $TARGETARCH" >&2; exit 1; \
    fi; \
    mosdepth --version

RUN python3 -m pip install --upgrade setuptools --break-system-packages
# RUN python3 -m pip install boto3

RUN git clone https://github.com/shinichinamba/ascairn.git && \
	cd ascairn && \
	python3 -m pip install . --break-system-packages && \
	which ascairn 

