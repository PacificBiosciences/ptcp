#!/bin/bash

set -euo pipefail

MYDIR=$(readlink -f "$(dirname "$0")")
source "${MYDIR}/build_versions.sh"

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  --help        Show this help message"
    echo "  --artifactory Build image with Artifactory registry prefix"
    echo "  --dryrun      Print build commands without executing them"
    echo "  --save-tarball Save Docker image as a compressed tarball"
    exit 1
}

if ! command -v docker &> /dev/null; then
    echo "Error: docker is not installed or not in PATH"
    exit 1
 fi

# Image settings
IMAGE_TAG=latest
IMAGE_FULL_TAG="${IMAGE_NAME}:${IMAGE_TAG}"

# Artifacory settings
ARTIFACTORY_URL=artifactory.pacificbiosciences.com
DOCKER_REPO=docker
ORG_NAME=cb

DRYRUN=false
SAVE_TARBALL=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)
            usage
            ;;
        --artifactory)
            IMAGE_FULL_TAG="${ARTIFACTORY_URL}/${DOCKER_REPO}/${ORG_NAME}/${IMAGE_NAME}:${IMAGE_TAG}"
            shift
            ;;
        --dryrun)
            DRYRUN=true
            shift
            ;;
        --save-tarball)
            SAVE_TARBALL=true
            shift
            ;;
        *)
            echo "Error: Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

echo "Building ${IMAGE_FULL_TAG}"

BUILD_CMD="docker buildx build \
  --rm \
  --platform linux/amd64 \
  --build-arg RUST_VERSION=${RUST_VERSION} \
  --build-arg SMRTLINK_VERSION=${SMRTLINK_VERSION} \
  --build-arg PARAPHASE_COMMIT=${PARAPHASE_COMMIT} \
  --build-arg TRGT_VERSION=${TRGT_VERSION} \
  --build-arg SAWFISH_VERSION=${SAWFISH_VERSION} \
  --build-arg PTCPQC_VERSION=${PTCPQC_VERSION} \
  -t \"${IMAGE_FULL_TAG}\" \
  -f \"${MYDIR}/Dockerfile\" \
  ."

if [ "$DRYRUN" = true ]; then
    echo "Dry run - would execute:"
    echo "$BUILD_CMD"
else
    eval "$BUILD_CMD"
    echo "Build completed successfully!"

    if [ "$SAVE_TARBALL" = true ]; then
        echo "Saving Docker tarball..."
        docker save "${IMAGE_FULL_TAG}" | gzip -c > "${IMAGE_NAME}_${IMAGE_TAG}.tar.gz"
        echo "Saved tarball as ${IMAGE_NAME}_${IMAGE_TAG}.tar.gz"
    fi
fi
