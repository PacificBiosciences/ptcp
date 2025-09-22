FROM ubuntu:focal

LABEL maintainer="Billy Rowell <wrowell@pacb.com>, Tom Mokveld <tmokveld@pacificbiosciences>, Guilherme De Sena Brandine <gbrandine@pacificbiosciences.com>"

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
  git \
  curl \
  python3 \
  python3-pip \
  && rm -rf /var/lib/apt/lists/*

ARG RUST_VERSION
ENV RUSTUP_HOME=/usr/local/rustup \
  CARGO_HOME=/usr/local/cargo \
  PATH=/usr/local/cargo/bin:$PATH \
  RUST_VERSION=${RUST_VERSION}
RUN curl https://sh.rustup.rs -sSf | sh -s -- -q -y --default-toolchain=$RUST_VERSION && \
  chmod -R a+w $RUSTUP_HOME $CARGO_HOME

RUN locale-gen en_US.UTF-8
RUN mkdir -p /opt/pacbio
RUN mkdir -p /usr/local/bin/

ARG SMRTLINK_VERSION
RUN wget --progress=dot:mega \
  https://downloads.pacbcloud.com/public/software/installers/smrtlink-release_${SMRTLINK_VERSION}.zip \
  && unzip smrtlink-release_${SMRTLINK_VERSION}.zip \
  && INSTALLER=$(find . -name "smrtlink-release_*.run" | head -1) \
  && ${INSTALLER} --batch --rootdir /opt/pacbio/smrttools --smrttools-only \
  && rm -f smrtlink-release_${SMRTLINK_VERSION}.zip ${INSTALLER}
ENV SMRT_ROOT "/opt/pacbio/smrttools"
ENV PATH "/opt/pacbio/smrttools/smrtcmds/bin:/opt/pacbio/smrttools/smrtcmds/developer/bin:/opt/pacbio/smrttools/current/bundles/smrttools/current/private/otherbins/all/bin:${PATH}"

ARG TRGT_VERSION
RUN mkdir -p /tmp/trgt_install \
  && cd /tmp/trgt_install \
  && wget --progress=dot:mega \
  https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
  && tar -xzf trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
  && cp trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /opt/pacbio/smrttools/install/smrtlink-release_${SMRTLINK_VERSION}/bundles/smrttools/install/smrttools-release_${SMRTLINK_VERSION}/private/pacbio/trgt/bin/trgt \
  && cd / \
  && rm -rf /tmp/trgt_install

ARG PARAPHASE_COMMIT
RUN git clone https://github.com/PacificBiosciences/paraphase.git /opt/paraphase \
  && cd /opt/paraphase \
  && git checkout ${PARAPHASE_COMMIT} \
  && python3 setup.py install \
  && pip3 install -e . \
  && cd /

# NOTE: This is temporary, for now ptcp will ship along a binary of the latest version of ptcp-qc 
ARG PTCPQC_ARCHIVE
COPY binaries/${PTCPQC_ARCHIVE} /tmp/
RUN mkdir -p /tmp/ptcp-qc-extract && \
  tar -xzf /tmp/${PTCPQC_ARCHIVE} -C /tmp/ptcp-qc-extract && \
  find /tmp/ptcp-qc-extract -type f -executable -exec cp {} /usr/local/bin/ptcp-qc \; && \
  chmod +x /usr/local/bin/ptcp-qc && \
  rm -rf /tmp/ptcp-qc-extract /tmp/${PTCPQC_ARCHIVE}

ARG SAWFISH_VERSION
RUN mkdir -p /tmp/sawfish_install \
  && cd /tmp/sawfish_install \
  && wget --progress=dot:mega \
  https://github.com/PacificBiosciences/sawfish/releases/download/${SAWFISH_VERSION}/sawfish-${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
  && tar -xzf sawfish-${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
  && cp sawfish-${SAWFISH_VERSION}-x86_64-unknown-linux-gnu/bin/sawfish /usr/local/bin/ \
  && cd / \
  && rm -rf /tmp/sawfish_install

COPY scripts/* /opt/scripts/
COPY scripts/havanno/havanno /opt/scripts/havanno
COPY scripts/smn-homology/smn_homology /opt/scripts/smn_homology
RUN chmod +x /opt/scripts/*

RUN groupadd -g 1000 bfxuser && \
  useradd -u 1000 -g 1000 -m -s /bin/bash bfxuser

ENV PATH="${PATH}:/opt/scripts"
RUN echo "PATH=${PATH}" >> /etc/bash.bashrc
