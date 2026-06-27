# Stage 1: Compile the Rust extension and produce a wheel.
FROM python:3.12 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

RUN pip install --no-cache-dir maturin

WORKDIR /build/ribo_rs
COPY ribo_rs/ .
# maturin build --release produces a platform wheel in target/wheels/.
# (build_rust.sh uses `maturin develop` which installs into the active venv directly
# and doesn't survive the stage boundary — `maturin build` is the right choice here.)
RUN maturin build --release


# Stage 2: Slim runtime image.
FROM python:3.12-slim

# Java JRE for VARNA headless renderer + font/harfbuzz libs it needs at runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        default-jre-headless \
        libharfbuzz0b \
        fontconfig \
        fonts-dejavu \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Rust extension wheel from builder stage.
COPY --from=builder /build/ribo_rs/target/wheels/ /tmp/wheels/
RUN pip install --no-cache-dir /tmp/wheels/ribo_rs-*.whl && rm -rf /tmp/wheels

# Python package (ribo_switch) + service and viz dependencies.
# Deliberately excludes ribo_rs/ source so it never shadows the installed wheel above.
COPY pyproject.toml README.md ./
COPY ribo_switch/ ribo_switch/
RUN pip install --no-cache-dir ".[service,viz]"

# VARNA headless renderer JAR.
COPY deps/VARNAv3-93.jar /app/deps/VARNAv3-93.jar
ENV VARNA_JAR=/app/deps/VARNAv3-93.jar

EXPOSE 8001
CMD ["uvicorn", "ribo_switch.service:app", "--host", "0.0.0.0", "--port", "8001"]
