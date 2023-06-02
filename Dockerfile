FROM ubuntu:20.04 AS builder
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase tabix \
     g++ cmake make automake autoconf bash-completion pigz zlib1g-dev

ADD https://github.com/ebiggers/libdeflate/archive/v1.18.tar.gz /tmp
RUN tar xzf /tmp/v1.18.tar.gz -C /tmp
WORKDIR /tmp/libdeflate-1.18
RUN cmake -B build && cmake --build build -j$(nproc) && cmake --install build && cp /usr/local/lib/libdeflate.so.0 /lib/x86_64-linux-gnu/libdeflate.so.0

ADD . /src
WORKDIR /src
RUN rm -f CMakeCache.txt && cmake -DCMAKE_BUILD_TYPE=Release /src && make clean && make -j$(nproc)
RUN ctest -V


FROM ubuntu:20.04
RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     bcftools less
COPY --from=builder /src/spvcf /usr/local/bin/spvcf
COPY --from=builder /src/external/src/htslib/bgzip /usr/local/bin/bgzip
COPY --from=builder /src/external/src/htslib/tabix /usr/local/bin/tabix
COPY --from=builder /usr/local/lib/libdeflate.so.0 /lib/x86_64-linux-gnu/libdeflate.so.0
CMD ["spvcf"]
