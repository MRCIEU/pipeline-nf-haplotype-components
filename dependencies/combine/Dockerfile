FROM alpine:3.23 AS build

RUN apk add --update --no-cache \
        build-base \
        cmake \
        curl \
        git \
        openssh

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs/ | sh -s -- -y

RUN git clone https://github.com/MRCIEU/combine-pbwt-sparsepaint.git && \
    cd combine-pbwt-sparsepaint && \
    git checkout v0.3.0

WORKDIR /combine-pbwt-sparsepaint

RUN source $HOME/.cargo/env && \
    cargo build --release

FROM alpine:3.23

RUN apk add --update --no-cache \
        procps

COPY --from=build /combine-pbwt-sparsepaint/target/release/combine /usr/local/bin/combine

CMD ["combine", "-h"]
