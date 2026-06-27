# Stage 1: Compile the Rust extension and produce a wheel.
# Uses the full python:3.12 image (has gcc/dev headers) plus the Rust toolchain.
FROM python:3.12 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

RUN pip install --no-cache-dir maturin

# Copy only the Rust crate; avoids dragging the rest of the repo into this layer.
WORKDIR /build/ribo_rs
COPY ribo_rs/ .
# maturin build --release produces a platform wheel in target/wheels/.
# build_rust.sh uses `maturin develop --release` which installs into the active venv
# directly — that install doesn't survive the stage boundary. Using `maturin build`
# instead creates a portable wheel we can pip-install in the final stage.
RUN maturin build --release


# Stage 2: Slim runtime image.
FROM python:3.12-slim

WORKDIR /app

# Install the compiled Rust extension from the wheel produced above.
# Inter-stage COPY is unaffected by .dockerignore, so target/wheels/ is always fresh.
COPY --from=builder /build/ribo_rs/target/wheels/ /tmp/wheels/
RUN pip install --no-cache-dir /tmp/wheels/ribo_rs-*.whl && rm -rf /tmp/wheels

# Copy only the Python package source — deliberately excludes ribo_rs/ so the source
# directory never shadows the site-packages install of ribo_rs above.
COPY pyproject.toml README.md ./
COPY ribo_switch/ ribo_switch/
RUN pip install --no-cache-dir ".[service]"

EXPOSE 8001
CMD ["uvicorn", "ribo_switch.service:app", "--host", "0.0.0.0", "--port", "8001"]
