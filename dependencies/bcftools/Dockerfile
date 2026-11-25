FROM gcc:15 AS build

RUN wget https://github.com/samtools/bcftools/releases/download/1.24/bcftools-1.24.tar.bz2 && \
    tar -xf bcftools-1.24.tar.bz2 && \
    cd bcftools-1.24 && \
    make && \
    cd htslib-1.24 && \
    make

FROM debian:13-slim

RUN apt-get -y update && \
        apt-get -y install \
        curl \
        procps

COPY --from=build /bcftools-1.24/bcftools usr/local/bin/bcftools
COPY --from=build /bcftools-1.24/htslib-1.24/tabix usr/local/bin/tabix
COPY --from=build /bcftools-1.24/htslib-1.24/bgzip usr/local/bin/bgzip
COPY --from=build /bcftools-1.24/LICENSE /usr/share/licenses/bcftools/LICENSE
COPY --from=build /bcftools-1.24/htslib-1.24/LICENSE /usr/share/licenses/htslib/LICENSE

CMD ["bcftools", "-h"]
