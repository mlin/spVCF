FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && \
     apt-get -qq install -y --no-install-recommends --no-install-suggests \
     curl wget ca-certificates git-core less netbase \
     g++ cmake make pigz

ADD . /src
WORKDIR /src
RUN rm -f CMakeCache.txt && cmake -DCMAKE_BUILD_TYPE=Release /src && make clean && make -j$(nproc)
CMD ctest -V
