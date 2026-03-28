# RuvLTRA Training Pipeline
# Supports: LoRA SFT, DPO, imatrix calibration, GGUF conversion
# Target: Cloud Run Jobs with L4 GPU or Vertex AI

FROM nvidia/cuda:12.4.1-devel-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-venv python3.11-dev python3-pip \
    git cmake build-essential curl \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

# Build llama.cpp with CUDA support for imatrix + GGUF conversion
RUN git clone --depth 1 https://github.com/ggerganov/llama.cpp /opt/llama.cpp \
    && cd /opt/llama.cpp \
    && cmake -B build -DGGML_CUDA=ON -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build --config Release -j$(nproc) --target llama-imatrix llama-quantize

# --- Runtime stage ---
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PATH="/opt/llama.cpp/build/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-venv python3-pip git \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

COPY --from=builder /opt/llama.cpp/build/bin/llama-imatrix /opt/llama.cpp/build/bin/llama-imatrix
COPY --from=builder /opt/llama.cpp/build/bin/llama-quantize /opt/llama.cpp/build/bin/llama-quantize

RUN pip install --no-cache-dir \
    torch==2.3.1 \
    transformers>=4.44.0 \
    peft>=0.12.0 \
    trl>=0.9.0 \
    datasets>=2.20.0 \
    huggingface_hub>=0.24.0 \
    llama-cpp-python>=0.2.80 \
    accelerate>=0.33.0 \
    bitsandbytes>=0.43.0 \
    sentencepiece \
    protobuf \
    safetensors

WORKDIR /app
COPY scripts/training/ /app/

ENTRYPOINT ["python", "-u"]
CMD ["run_calibration.py"]
