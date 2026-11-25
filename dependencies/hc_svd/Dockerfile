FROM rust:1.94-bookworm AS build

RUN git clone https://github.com/MRCIEU/hc-svd.git && \
    cd hc-svd && \
    git checkout v0.2.0

WORKDIR /hc-svd

RUN cargo build --release --locked

# Debian slim (rather than scratch or Alpine) final stage: Nextflow needs a shell and ps in the container,
# and glibc's allocator is faster than musl's for this kind of work
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*

COPY --from=build /hc-svd/target/release/hcsvd /usr/local/bin/hcsvd
COPY --from=build /hc-svd/LICENSE /hc-svd/NOTICE /usr/share/licenses/hc-svd/

CMD ["hcsvd", "--help"]
