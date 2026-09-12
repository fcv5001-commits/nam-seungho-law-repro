"""Open-data download helpers for Nam Seung-ho Law reproducibility checks."""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path


SOURCES = {
    "pantheon_plus": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat",
    "planck_tt_binned": "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-binned_R3.01.txt",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(name: str, cache_dir: str | Path = "data") -> tuple[Path, str]:
    """Download a named public dataset once and return path and SHA-256."""
    if name not in SOURCES:
        raise KeyError(f"unknown source: {name}")
    directory = Path(cache_dir)
    directory.mkdir(parents=True, exist_ok=True)
    suffix = ".dat" if name == "pantheon_plus" else ".txt"
    target = directory / f"{name}{suffix}"
    if not target.exists():
        request = urllib.request.Request(
            SOURCES[name], headers={"User-Agent": "nam-seungho-law-repro/1.0"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            target.write_bytes(response.read())
    return target, sha256(target)


def write_manifest(cache_dir: str | Path = "data") -> Path:
    directory = Path(cache_dir)
    entries = {}
    for name, url in SOURCES.items():
        path, checksum = fetch(name, directory)
        entries[name] = {"url": url, "file": str(path), "sha256": checksum}
    output = directory / "manifest.json"
    output.write_text(json.dumps(entries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return output


if __name__ == "__main__":
    print(write_manifest())
