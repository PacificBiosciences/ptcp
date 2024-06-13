FROM ubuntu:focal

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
  && apt-get -qq install \
    wget \
    unzip \
    rsync \
    locales \
    build-essential \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    jq \
  && rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US.UTF-8
RUN mkdir -p /opt/pacbio
RUN mkdir -p /usr/local/bin/

ARG SMRTLINK_VERSION
RUN wget --progress=dot:mega \
  https://downloads.pacbcloud.com/public/software/installers/smrtlink-release_${SMRTLINK_VERSION}.zip
RUN unzip smrtlink-release_${SMRTLINK_VERSION}.zip
RUN INSTALLER=$(find . -name "smrtlink-release_*.run" | head -1) && \
    ${INSTALLER} --batch \
    --rootdir /opt/pacbio/smrttools \
    --smrttools-only
RUN rm -f smrtlink-release_${SMRTLINK_VERSION}.zip ${INSTALLER}
ENV SMRT_ROOT "/opt/pacbio/smrttools"
ENV PATH "/opt/pacbio/smrttools/smrtcmds/bin:/opt/pacbio/smrttools/smrtcmds/developer/bin:/opt/pacbio/smrttools/current/bundles/smrttools/current/private/otherbins/all/bin:${PATH}"

ARG TRGT_VERSION
RUN wget --progress=dot:mega \
  https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN tar --no-same-owner -zxvf  trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
  && rm trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN mv trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /opt/pacbio/smrttools/install/smrtlink-release_${SMRTLINK_VERSION}/bundles/smrttools/install/smrttools-release_${SMRTLINK_VERSION}/private/pacbio/trgt/bin/trgt \
  && chmod +x /opt/pacbio/smrttools/install/smrtlink-release_${SMRTLINK_VERSION}/bundles/smrttools/install/smrttools-release_${SMRTLINK_VERSION}/private/pacbio/trgt/bin/trgt
RUN rm /opt/pacbio/smrttools/install/smrtlink-release_${SMRTLINK_VERSION}/bundles/smrttools/install/smrttools-release_${SMRTLINK_VERSION}/private/pacbio/trgt/bin/trvz \
  && rm /opt/pacbio/smrttools/install/smrtlink-release_${SMRTLINK_VERSION}/bundles/smrttools/install/smrttools-release_${SMRTLINK_VERSION}/private/pacbio/trgt/binwrap/trvz

ARG BEDTOOLS_VERSION
RUN wget --progress=dot:mega \
  https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static \
  -O /usr/local/bin/bedtools \
  && chmod +x /usr/local/bin/bedtools

ARG MINIPILEUP_VERSION
RUN wget --progress=dot:mega \
  https://github.com/lh3/minipileup/archive/${MINIPILEUP_VERSION}.tar.gz \
  && tar --no-same-owner -zxvf ${MINIPILEUP_VERSION}.tar.gz --directory /opt \
  && rm ${MINIPILEUP_VERSION}.tar.gz
RUN cd /opt/minipileup-${MINIPILEUP_VERSION} \
  && make \
  && ln -s /opt/minipileup-${MINIPILEUP_VERSION}/minipileup /usr/local/bin/

COPY scripts/* /opt/scripts/
RUN chmod +x /opt/scripts/*

ENV PATH="${PATH}:/opt/scripts"
RUN echo "PATH=${PATH}" >> /etc/bash.bashrc