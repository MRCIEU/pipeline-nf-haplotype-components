FROM alpine:3.23 AS build

RUN apk update && \
    apk add \
        git \
        build-base \
        zlib-dev \
        bzip2-dev \
        xz-dev \
        libdeflate \
        curl-dev \
        autoconf \
        automake \
        libcrypto3

RUN git clone https://github.com/samtools/htslib && \
    cd htslib && \
    git checkout 1.23 && \
    git submodule update --init --recursive  && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

RUN git clone https://github.com/richarddurbin/pbwt && \
    cd pbwt && \
    git checkout v2.1 && \
    sed -i '5s/$/ -lcrypto/' Makefile && \
    make


FROM alpine:3.23

RUN apk update && \
    apk add \
        bzip2-dev \
        xz-dev \
        curl-dev \
        procps

COPY --from=build pbwt/pbwt usr/local/bin/pbwt

CMD ["pbwt"]
