#!/usr/bin/env python3
"""Generate a synthetic FASTQ with realistic entropy (so it compresses ~4x, like real data).

Record layout is fixed-width, 316 bytes:
  @r<9 digits>\n  <150 bp>\n  +\n  <150 qual>\n
Sequence is uniform random ACGT.  Quality uses the 4-bin NovaSeq alphabet with a
position-dependent profile that degrades toward the 3' end, which is what gives real
FASTQ its compressibility.
"""
import argparse, sys
import numpy as np

RL = 150
HDR = 12          # '@r' + 9 digits + '\n'
REC = HDR + RL + 1 + 2 + RL + 1   # 316

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
# NovaSeq-style binned quality scores: Q2 '#', Q12 '-', Q23 '8', Q37 'F'
QBINS = np.frombuffer(b"#-8F", dtype=np.uint8)


def qual_profile(rl):
    """P(bin) per cycle: high quality early, degrading late."""
    x = np.linspace(0.0, 1.0, rl)
    p_hi = 0.97 - 0.35 * x**2          # fraction in the top bin
    p = np.empty((rl, 4))
    p[:, 3] = p_hi
    p[:, 2] = (1 - p_hi) * 0.55
    p[:, 1] = (1 - p_hi) * 0.33
    p[:, 0] = (1 - p_hi) * 0.12
    return p / p.sum(axis=1, keepdims=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=250_000)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    prof = qual_profile(RL)
    # Thresholds on a 0-255 scale so the draw can come straight from raw RNG bytes,
    # which is far cheaper than drawing float64 uniforms.
    cdf = np.cumsum(prof, axis=1)[:, :3].T.copy()   # (3, RL)
    thr = np.clip(np.rint(cdf * 256), 0, 255).astype(np.uint8)

    pow10 = (10 ** np.arange(8, -1, -1)).astype(np.int64)

    buf = np.empty((args.chunk, REC), dtype=np.uint8)
    buf[:, 0] = ord("@")
    buf[:, 1] = ord("r")
    buf[:, HDR - 1] = ord("\n")
    buf[:, HDR + RL] = ord("\n")
    buf[:, HDR + RL + 1] = ord("+")
    buf[:, HDR + RL + 2] = ord("\n")
    buf[:, REC - 1] = ord("\n")

    written = 0
    with open(args.out, "wb", buffering=1 << 20) as fh:
        while written < args.reads:
            n = min(args.chunk, args.reads - written)
            v = buf[:n]
            idx = np.arange(written, written + n, dtype=np.int64)
            v[:, 2:HDR - 1] = ((idx[:, None] // pow10) % 10 + 48).astype(np.uint8)
            raw = np.frombuffer(rng.bytes(n * 2 * RL), dtype=np.uint8).reshape(n, 2 * RL)
            v[:, HDR:HDR + RL] = BASES[raw[:, :RL] & 3]
            u = raw[:, RL:]
            qi = (u >= thr[0]).astype(np.uint8) + (u >= thr[1]) + (u >= thr[2])
            v[:, HDR + RL + 3:REC - 1] = QBINS[qi]
            fh.write(v.tobytes())
            written += n
            if written % 25_000_000 == 0:
                print(f"  {written:,} reads", file=sys.stderr, flush=True)
    print(f"done: {written:,} reads, {written * REC:,} bytes", file=sys.stderr)


if __name__ == "__main__":
    main()
