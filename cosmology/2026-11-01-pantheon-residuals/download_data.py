"""Fetch official Pantheon+ inputs and require the recorded SHA256 checksums."""
from pathlib import Path
import hashlib,json,urllib.request

ROOT=Path(__file__).parent
BASE='https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/'
EXPECTED=json.loads((ROOT/'reference_results/summary.json').read_text())['input_sha256']

def main():
    directory=ROOT/'data';directory.mkdir(exist_ok=True)
    for name,digest in EXPECTED.items():
        target=directory/name
        if not target.exists():
            url=BASE+name.replace('+','%2B')
            with urllib.request.urlopen(url,timeout=120) as response:
                data=response.read()
            if hashlib.sha256(data).hexdigest()!=digest:
                raise RuntimeError(f'Official bytes differ from recorded input: {name}; no file installed')
            target.write_bytes(data)
        if hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
            raise RuntimeError(f'Local input checksum mismatch: {name}; do not continue')
        print(f'INPUT HASH MATCH {name}')

if __name__=='__main__':main()
