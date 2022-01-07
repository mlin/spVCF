FROM ubuntu:18.04 AS builder
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase tabix \
     g++ cmake make automake autoconf bash-completion pigz zlib1g-dev

ADD https://github.com/ebiggers/libdeflate/archive/v1.8.tar.gz /tmp
RUN tar xzf /tmp/v1.8.tar.gz -C /tmp
WORKDIR /tmp/libdeflate-1.8
RUN make -j $(nproc) && make install

ADD . /src
WORKDIR /src
RUN rm -f CMakeCache.txt && cmake -DCMAKE_BUILD_TYPE=Release /src && make clean && make -j$(nproc)
RUN ctest -V


FROM ubuntu:20.04
RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     tabix bcftools less
COPY --from=builder /src/spvcf /usr/local/bin/spvcf
CMD ["spvcf"]
