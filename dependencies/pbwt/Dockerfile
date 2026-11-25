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

# pbwt has no LICENSE file: its source files state the Apache License 2.0
RUN wget -q -O /pbwt-LICENSE https://www.apache.org/licenses/LICENSE-2.0.txt && \
    printf '%s\n' \
        'pbwt v2.1 (https://github.com/richarddurbin/pbwt)' \
        'Copyright (C) Genome Research Limited, 2013-' \
        'Licensed under the Apache License, Version 2.0 (see LICENSE), as stated in its source files.' \
        'It is statically linked with htslib 1.23 (see ../htslib/LICENSE).' > /pbwt-NOTICE


FROM alpine:3.23

RUN apk update && \
    apk add \
        bzip2-dev \
        xz-dev \
        curl-dev \
        procps

COPY --from=build pbwt/pbwt usr/local/bin/pbwt
COPY --from=build /pbwt-LICENSE /usr/share/licenses/pbwt/LICENSE
COPY --from=build /pbwt-NOTICE /usr/share/licenses/pbwt/NOTICE
COPY --from=build /htslib/LICENSE /usr/share/licenses/htslib/LICENSE

CMD ["pbwt"]
