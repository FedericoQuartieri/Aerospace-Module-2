# Development environment for the project: OpenFOAM-13 + Mutation++.
#
# Ubuntu 24.04 "noble": dl.openfoam.org publishes openfoam13 for arm64 too,
# so on Apple Silicon the image runs natively rather than under emulation.
# This matters: the run times are the real ones of the machine.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# ---- OpenFOAM-13 from the OpenFOAM Foundation repository
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates wget gnupg \
 && wget -qO /etc/apt/trusted.gpg.d/openfoam.asc https://dl.openfoam.org/gpg.key \
 && echo "deb http://dl.openfoam.org/ubuntu noble main" > /etc/apt/sources.list.d/openfoam.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends openfoam13 \
 && rm -rf /var/lib/apt/lists/*

# ---- tools: Mutation++ build and comparison scripts
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential cmake git patch \
        python3 python3-numpy python3-matplotlib \
        less nano time \
 && rm -rf /var/lib/apt/lists/*

# ---- Mutation++ inside the image, not on the bind mount
# the project would look for it in thirdParty/Mutationpp (see etc/config.sh/mutationpp):
# the entrypoint redirects MPP_DIRECTORY here, so the build does not go through
# the filesystem shared with macOS (slow) and the image is self-contained
#
# the revision is pinned for reproducibility: the image rebuilds identically
# months later. Mutation++ is installed as it comes from upstream, without
# patches: the project links against the library as it comes from upstream
#
# -j2 and not -j$(nproc): the units that include Eigen take more than 1 GB
# each and the Docker Desktop VM has 4 in total, so -j8 runs out of memory.
ENV MPP_DIRECTORY=/opt/Mutationpp
ENV MPP_REVISION=e8edf4f3a7f2bc22b4c62a23060dbe86967f1b3e

RUN git clone https://github.com/mutationpp/Mutationpp $MPP_DIRECTORY \
 && git -C $MPP_DIRECTORY checkout --detach $MPP_REVISION \
 && cmake -S $MPP_DIRECTORY -B $MPP_DIRECTORY/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_TESTING=OFF \
        -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER \
        -DCMAKE_INSTALL_PREFIX=$MPP_DIRECTORY/install \
 && cmake --build $MPP_DIRECTORY/build -j2 --target install \
 && rm -rf $MPP_DIRECTORY/build

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

WORKDIR /project
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
