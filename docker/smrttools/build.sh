#!/bin/bash

MYDIR=$(readlink -f "$(dirname "$0")")
SMRTLINK_VERSION=13.1.0.221970
TRGT_VERSION=1.0.0
BEDTOOLS_VERSION=2.31.0
MINIPILEUP_VERSION=78d58cfeed6ff32b721af19d6ede79f3b07cfe61
# Image info
IMAGE_NAME=smrttools
IMAGE_TAG=${SMRTLINK_VERSION}DIRTY

docker build \
  --rm \
  --build-arg SMRTLINK_VERSION=${SMRTLINK_VERSION} \
  --build-arg TRGT_VERSION=${TRGT_VERSION} \
  --build-arg BEDTOOLS_VERSION=${BEDTOOLS_VERSION} \
  --build-arg MINIPILEUP_VERSION=${MINIPILEUP_VERSION} \
  --build-arg IMAGE_NAME=${IMAGE_NAME} \
  --build-arg IMAGE_TAG=${IMAGE_TAG} \
  -t "smrttools:${IMAGE_TAG}" \
  -f ${MYDIR}/Dockerfile \
  $(pwd)

# sha256sum=$(docker inspect --format='{{index .RepoDigests 0}}' "smrttools:${IMAGE_TAG}" 2>/dev/null | sha256sum | cut -d ' ' -f 1)

# to save the docker tarball for DNAnexus
# docker save "smrttools:${IMAGE_TAG}" | gzip -c > "smrttools_${IMAGE_TAG}.tar.gz"

# to build the sif image for singularity, internal testing
# singularity build "docker___quay.io_pacbio_smrttools@sha256_${sha256sum}.sif" docker-daemon://smrttools:${IMAGE_TAG}
