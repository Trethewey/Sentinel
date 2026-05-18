FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools bcftools \
        build-essential libssl-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/sentinel
COPY pyproject.toml README.md ./
COPY src ./src

RUN pip install --no-cache-dir .

ENTRYPOINT ["sentinel"]
CMD ["--help"]
