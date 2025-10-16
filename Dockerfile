FROM ghcr.io/astral-sh/uv:python3.12-trixie-slim

LABEL name="ms2rescore"

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
    && useradd --system --gid 999 --uid 999 --create-home nonroot

# Install the project into `/ms2rescore`
WORKDIR /ms2rescore

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Ensure installed tools can be executed out of the box
ENV UV_TOOL_BIN_DIR=/usr/local/bin

RUN apt-get update && apt-get install -y procps

# Then, add the rest of the project source code and install it
# Installing separately from its dependencies allows optimal layer caching
ADD pyproject.toml /ms2rescore/pyproject.toml
ADD LICENSE /ms2rescore/LICENSE
ADD README.md /ms2rescore/README.md
ADD MANIFEST.in /ms2rescore/MANIFEST.in
ADD uv.lock /ms2rescore/uv.lock
ADD ms2rescore /ms2rescore/ms2rescore

# Install the project and its dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --no-dev

# Place executables in the environment at the front of the path
ENV PATH="/ms2rescore/.venv/bin:$PATH"

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []

# Use the non-root user to run our application
USER nonroot
