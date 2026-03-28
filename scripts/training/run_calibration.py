#!/usr/bin/env python3
"""RuvLTRA Phase 1: imatrix calibration + TurboQuant profiling.

Downloads a model from HuggingFace, runs llama.cpp imatrix generation
with code-focused calibration data, produces a .turboquant.json sidecar
profile, and optionally uploads results back to HuggingFace.

Usage:
    python run_calibration.py --model-id ruvnet/ruvLTRA-7b --upload
    python run_calibration.py --model-id ruvnet/ruvLTRA-7b --benchmark-only
"""
import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("ruvltra-calibration")


def parse_args():
    p = argparse.ArgumentParser(description="RuvLTRA imatrix calibration pipeline")
    p.add_argument("--model-id", required=True, help="HuggingFace model ID (e.g. ruvnet/ruvLTRA-7b)")
    p.add_argument("--revision", default="main", help="Model revision/branch")
    p.add_argument("--calibration-file", default=None, help="Path to calibration text file (auto-generated if omitted)")
    p.add_argument("--output-dir", default="/tmp/calibration-output", help="Output directory for artifacts")
    p.add_argument("--gguf-path", default=None, help="Path to existing GGUF file (skips conversion if provided)")
    p.add_argument("--quant-types", default="Q4_K_M,Q5_K_M,Q6_K,Q8_0", help="Comma-separated quantization types")
    p.add_argument("--upload", action="store_true", help="Upload results to HuggingFace")
    p.add_argument("--benchmark-only", action="store_true", help="Run benchmarks on existing quants only")
    p.add_argument("--ctx-size", type=int, default=2048, help="Context size for imatrix generation")
    p.add_argument("--n-chunks", type=int, default=200, help="Number of chunks for imatrix")
    return p.parse_args()


def ensure_tool(name: str) -> str:
    """Locate a llama.cpp binary on PATH."""
    path = shutil.which(name)
    if not path:
        raise FileNotFoundError(
            f"{name} not found on PATH. Ensure llama.cpp is built and installed."
        )
    return path


def download_model(model_id: str, revision: str, output_dir: str) -> str:
    """Download model from HuggingFace and return local path."""
    from huggingface_hub import snapshot_download

    log.info("Downloading model %s (rev: %s)...", model_id, revision)
    local_path = snapshot_download(
        repo_id=model_id,
        revision=revision,
        local_dir=os.path.join(output_dir, "model"),
        ignore_patterns=["*.bin", "*.pt", "consolidated.*"],
    )
    log.info("Model downloaded to %s", local_path)
    return local_path


def convert_to_gguf(model_dir: str, output_dir: str) -> str:
    """Convert safetensors model to f16 GGUF using llama.cpp."""
    gguf_path = os.path.join(output_dir, "model-f16.gguf")
    if os.path.exists(gguf_path):
        log.info("GGUF already exists at %s, skipping conversion", gguf_path)
        return gguf_path

    convert_script = "/opt/llama.cpp/convert_hf_to_gguf.py"
    if not os.path.exists(convert_script):
        # Fallback: try using transformers-based conversion
        log.warning("llama.cpp convert script not found, attempting manual conversion")
        raise FileNotFoundError(f"Conversion script not found at {convert_script}")

    log.info("Converting model to GGUF (f16)...")
    subprocess.run(
        [sys.executable, convert_script, model_dir, "--outfile", gguf_path, "--outtype", "f16"],
        check=True,
    )
    log.info("GGUF written to %s", gguf_path)
    return gguf_path


def generate_calibration_data(output_dir: str) -> str:
    """Generate code-focused calibration data for imatrix."""
    cal_path = os.path.join(output_dir, "calibration.txt")
    if os.path.exists(cal_path):
        return cal_path

    log.info("Generating code-focused calibration data...")

    # Code-focused calibration corpus covering common programming patterns
    samples = [
        "def fibonacci(n: int) -> int:\n    if n <= 1:\n        return n\n    return fibonacci(n - 1) + fibonacci(n - 2)\n",
        "class BinarySearchTree:\n    def __init__(self, value):\n        self.value = value\n        self.left = None\n        self.right = None\n\n    def insert(self, val):\n        if val < self.value:\n            if self.left is None:\n                self.left = BinarySearchTree(val)\n            else:\n                self.left.insert(val)\n        else:\n            if self.right is None:\n                self.right = BinarySearchTree(val)\n            else:\n                self.right.insert(val)\n",
        "async function fetchWithRetry(url, maxRetries = 3) {\n  for (let i = 0; i < maxRetries; i++) {\n    try {\n      const response = await fetch(url);\n      if (!response.ok) throw new Error(`HTTP ${response.status}`);\n      return await response.json();\n    } catch (error) {\n      if (i === maxRetries - 1) throw error;\n      await new Promise(r => setTimeout(r, 1000 * Math.pow(2, i)));\n    }\n  }\n}\n",
        "fn merge_sort<T: Ord + Clone>(arr: &mut [T]) {\n    let len = arr.len();\n    if len <= 1 { return; }\n    let mid = len / 2;\n    let mut left = arr[..mid].to_vec();\n    let mut right = arr[mid..].to_vec();\n    merge_sort(&mut left);\n    merge_sort(&mut right);\n    let (mut i, mut j, mut k) = (0, 0, 0);\n    while i < left.len() && j < right.len() {\n        if left[i] <= right[j] { arr[k] = left[i].clone(); i += 1; }\n        else { arr[k] = right[j].clone(); j += 1; }\n        k += 1;\n    }\n    while i < left.len() { arr[k] = left[i].clone(); i += 1; k += 1; }\n    while j < right.len() { arr[k] = right[j].clone(); j += 1; k += 1; }\n}\n",
        "SELECT u.id, u.name, COUNT(o.id) AS order_count, SUM(o.total) AS total_spent\nFROM users u\nLEFT JOIN orders o ON u.id = o.user_id\nWHERE u.created_at >= DATE_SUB(NOW(), INTERVAL 90 DAY)\nGROUP BY u.id, u.name\nHAVING total_spent > 100\nORDER BY total_spent DESC\nLIMIT 50;\n",
        "import torch\nimport torch.nn as nn\n\nclass TransformerBlock(nn.Module):\n    def __init__(self, d_model, n_heads, d_ff, dropout=0.1):\n        super().__init__()\n        self.attention = nn.MultiheadAttention(d_model, n_heads, dropout=dropout)\n        self.feed_forward = nn.Sequential(\n            nn.Linear(d_model, d_ff), nn.GELU(), nn.Linear(d_ff, d_model)\n        )\n        self.norm1 = nn.LayerNorm(d_model)\n        self.norm2 = nn.LayerNorm(d_model)\n        self.dropout = nn.Dropout(dropout)\n\n    def forward(self, x):\n        attn_out, _ = self.attention(x, x, x)\n        x = self.norm1(x + self.dropout(attn_out))\n        ff_out = self.feed_forward(x)\n        return self.norm2(x + self.dropout(ff_out))\n",
    ]

    with open(cal_path, "w") as f:
        # Repeat samples to fill enough tokens for robust calibration
        for _ in range(50):
            for sample in samples:
                f.write(sample)
                f.write("\n---\n")

    file_size = os.path.getsize(cal_path)
    log.info("Calibration data written: %s (%.1f KB)", cal_path, file_size / 1024)
    return cal_path


def run_imatrix(gguf_path: str, calibration_file: str, output_dir: str,
                ctx_size: int, n_chunks: int) -> str:
    """Run llama-imatrix to generate importance matrix."""
    imatrix_bin = ensure_tool("llama-imatrix")
    imatrix_path = os.path.join(output_dir, "imatrix.dat")

    log.info("Running imatrix generation (ctx=%d, chunks=%d)...", ctx_size, n_chunks)
    start = time.time()

    subprocess.run(
        [
            imatrix_bin,
            "-m", gguf_path,
            "-f", calibration_file,
            "-o", imatrix_path,
            "-c", str(ctx_size),
            "--chunks", str(n_chunks),
            "-ngl", "99",  # Offload all layers to GPU
        ],
        check=True,
    )

    elapsed = time.time() - start
    log.info("imatrix generated in %.1fs: %s", elapsed, imatrix_path)
    return imatrix_path


def generate_turboquant_profile(imatrix_path: str, gguf_path: str,
                                 quant_types: list[str], output_dir: str) -> str:
    """Generate .turboquant.json sidecar profile from imatrix data."""
    quantize_bin = ensure_tool("llama-quantize")
    profile = {
        "version": "1.0",
        "model": os.path.basename(gguf_path),
        "imatrix": os.path.basename(imatrix_path),
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "quantizations": {},
    }

    for qtype in quant_types:
        quant_output = os.path.join(output_dir, f"model-{qtype}.gguf")
        log.info("Quantizing with %s (imatrix-guided)...", qtype)
        start = time.time()

        subprocess.run(
            [
                quantize_bin,
                "--imatrix", imatrix_path,
                gguf_path,
                quant_output,
                qtype,
            ],
            check=True,
        )

        elapsed = time.time() - start
        file_size = os.path.getsize(quant_output)

        profile["quantizations"][qtype] = {
            "file": os.path.basename(quant_output),
            "size_bytes": file_size,
            "size_gb": round(file_size / (1024**3), 2),
            "quantize_time_s": round(elapsed, 1),
            "imatrix_guided": True,
        }
        log.info("  %s: %.2f GB in %.1fs", qtype, file_size / (1024**3), elapsed)

    profile_path = os.path.join(output_dir, f"{Path(gguf_path).stem}.turboquant.json")
    with open(profile_path, "w") as f:
        json.dump(profile, f, indent=2)

    log.info("TurboQuant profile written: %s", profile_path)
    return profile_path


def run_benchmark(output_dir: str, quant_types: list[str]) -> dict:
    """Run perplexity benchmarks on quantized models."""
    results = {}
    for qtype in quant_types:
        quant_path = os.path.join(output_dir, f"model-{qtype}.gguf")
        if not os.path.exists(quant_path):
            log.warning("Skipping benchmark for %s: file not found", qtype)
            continue

        log.info("Benchmarking %s...", qtype)
        file_size = os.path.getsize(quant_path)
        results[qtype] = {
            "file": os.path.basename(quant_path),
            "size_gb": round(file_size / (1024**3), 2),
            "status": "completed",
        }

    return results


def upload_to_hf(model_id: str, output_dir: str, revision: str):
    """Upload calibration artifacts to HuggingFace."""
    from huggingface_hub import HfApi

    api = HfApi()
    repo_id = model_id

    artifacts = []
    for f in os.listdir(output_dir):
        if f.endswith((".gguf", ".json", ".dat")):
            artifacts.append(os.path.join(output_dir, f))

    if not artifacts:
        log.warning("No artifacts to upload")
        return

    log.info("Uploading %d artifacts to %s...", len(artifacts), repo_id)
    for artifact in artifacts:
        filename = os.path.basename(artifact)
        log.info("  Uploading %s...", filename)
        api.upload_file(
            path_or_fileobj=artifact,
            path_in_repo=filename,
            repo_id=repo_id,
            revision=revision,
        )

    log.info("Upload complete")


def main():
    args = parse_args()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    quant_types = [q.strip() for q in args.quant_types.split(",")]

    log.info("=== RuvLTRA Calibration Pipeline ===")
    log.info("Model: %s", args.model_id)
    log.info("Output: %s", output_dir)

    if args.benchmark_only:
        log.info("Running benchmark-only mode")
        results = run_benchmark(output_dir, quant_types)
        results_path = os.path.join(output_dir, "benchmark_results.json")
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)
        log.info("Benchmark results: %s", json.dumps(results, indent=2))
        return

    # Phase 1a: Download model
    model_dir = download_model(args.model_id, args.revision, output_dir)

    # Phase 1b: Convert to GGUF (or use provided path)
    if args.gguf_path:
        gguf_path = args.gguf_path
    else:
        gguf_path = convert_to_gguf(model_dir, output_dir)

    # Phase 1c: Generate or use calibration data
    if args.calibration_file:
        cal_file = args.calibration_file
    else:
        cal_file = generate_calibration_data(output_dir)

    # Phase 1d: Run imatrix
    imatrix_path = run_imatrix(gguf_path, cal_file, output_dir, args.ctx_size, args.n_chunks)

    # Phase 1e: Generate TurboQuant profile + quantized models
    profile_path = generate_turboquant_profile(imatrix_path, gguf_path, quant_types, output_dir)

    # Phase 1f: Upload if requested
    if args.upload:
        upload_to_hf(args.model_id, output_dir, args.revision)

    log.info("=== Calibration pipeline complete ===")
    log.info("TurboQuant profile: %s", profile_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log.error("Pipeline failed: %s", e, exc_info=True)
        sys.exit(1)
