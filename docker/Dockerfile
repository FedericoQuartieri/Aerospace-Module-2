# Ambiente di sviluppo del progetto: OpenFOAM-13 + Mutation++.
#
# Ubuntu 24.04 "noble": dl.openfoam.org pubblica openfoam13 anche per arm64,
# quindi su Apple Silicon l'immagine gira in nativo e non sotto emulazione.
# Questo conta: i tempi di calcolo sono quelli veri della macchina.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# ---- OpenFOAM-13 dal repository della OpenFOAM Foundation
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates wget gnupg \
 && wget -qO /etc/apt/trusted.gpg.d/openfoam.asc https://dl.openfoam.org/gpg.key \
 && echo "deb http://dl.openfoam.org/ubuntu noble main" > /etc/apt/sources.list.d/openfoam.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends openfoam13 \
 && rm -rf /var/lib/apt/lists/*

# ---- strumenti: build di Mutation++ e script di confronto
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential cmake git patch \
        python3 python3-numpy python3-matplotlib \
        less nano time \
 && rm -rf /var/lib/apt/lists/*

# ---- Mutation++ dentro l'immagine, non sul bind mount
# il progetto lo cercherebbe in thirdParty/Mutationpp (vedi etc/config.sh/mutationpp):
# l'entrypoint reindirizza MPP_DIRECTORY qui, cosi' la build non passa dal
# filesystem condiviso con macOS (lento) e l'immagine e' autosufficiente
#
# la revisione e' bloccata per riproducibilita': l'immagine si ricostruisce
# identica a distanza di mesi. Mutation++ e' installato come viene da monte,
# senza patch: il progetto si linka alla libreria come viene da monte
#
# -j2 e non -j$(nproc): le unita' che includono Eigen prendono oltre 1 GB
# ciascuna e la VM di Docker Desktop ne ha 4 in tutto, con -j8 si esaurisce.
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
