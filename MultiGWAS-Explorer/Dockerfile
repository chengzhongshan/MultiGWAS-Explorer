FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    OPEN_RESULT=0 \
    PIPELINE_ROOT=/opt/MultiGWAS-Explorer

WORKDIR /opt/MultiGWAS-Explorer

COPY . /opt/MultiGWAS-Explorer

RUN bash install/install_ubuntu.sh \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/MultiGWAS-Explorer/local/bin:/opt/MultiGWAS-Explorer/.venv-pipeline/bin:/opt/MultiGWAS-Explorer/local/perl5/bin:${PATH}"

CMD ["bash"]
