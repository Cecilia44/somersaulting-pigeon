#!/usr/bin/env python3
"""Download the official SAM 2.1 Hiera Large checkpoint."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlopen


URL = (
    "https://dl.fbaipublicfiles.com/segment_anything_2/092824/"
    "sam2.1_hiera_large.pt"
)
EXPECTED_BYTES = 898_083_611


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", default="checkpoints/sam2.1_hiera_large.pt", type=Path
    )
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_file() and output.stat().st_size == EXPECTED_BYTES:
        print(f"Checkpoint already exists: {output}")
        return 0

    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".part")
    print(f"Downloading {URL} ({EXPECTED_BYTES / 1024**2:.1f} MiB)")
    with urlopen(URL) as response, partial.open("wb") as handle:
        downloaded = 0
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)
            downloaded += len(chunk)
            print(f"\r{downloaded / EXPECTED_BYTES:6.1%}", end="", flush=True)
    print()
    if partial.stat().st_size != EXPECTED_BYTES:
        raise RuntimeError(
            f"Unexpected checkpoint size: {partial.stat().st_size} bytes; expected {EXPECTED_BYTES}"
        )
    partial.replace(output)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

