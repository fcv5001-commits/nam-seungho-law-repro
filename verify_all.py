"""Run every public reproduction script and hash the produced evidence."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


def run(script: str):
    subprocess.run([sys.executable, script], check=True)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    run("fit_pantheon_ns.py")
    run("cmb_residual_ns.py")
    outputs = sorted(Path("results").glob("*.json"))
    manifest = {str(path): digest(path) for path in outputs if path.name != "SHA256SUMS.json"}
    target = Path("results/SHA256SUMS.json")
    target.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
