FROM rust:1.94-alpine3.23 AS build

RUN apk add --update --no-cache \
        build-base \
        cmake \
        git

RUN git clone https://github.com/MRCIEU/combine-pbwt-sparsepaint.git && \
    cd combine-pbwt-sparsepaint && \
    git checkout v0.7.0

WORKDIR /combine-pbwt-sparsepaint

RUN cargo build --release --locked

# alpine (rather than scratch) final stage, because Nextflow needs a shell and ps in the container
FROM alpine:3.23

RUN apk add --update --no-cache \
        procps

COPY --from=build /combine-pbwt-sparsepaint/target/release/combine /usr/local/bin/combine
COPY --from=build /combine-pbwt-sparsepaint/LICENSE /combine-pbwt-sparsepaint/NOTICE /usr/share/licenses/combine/

CMD ["combine", "-h"]
