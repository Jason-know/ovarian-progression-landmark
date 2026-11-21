"""Restore the exact public input files used by the revision."""
from pathlib import Path
from urllib.request import Request, urlopen
import csv
import hashlib
import shutil

HERE = Path(__file__).resolve().parent


def digest(path, kind):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, kind).hexdigest()


def fetch():
    data = HERE / "data"
    data.mkdir(exist_ok=True)
    with (HERE / "results/source_manifest.csv").open(encoding="utf-8", newline="") as table:
        for row in csv.DictReader(table):
            path = data / row["file"]
            if path.exists() and digest(path, row["algorithm"]) == row["checksum"]:
                print(f"Verified {path.name}")
                continue
            temp = path.with_suffix(path.suffix + ".download")
            try:
                request = Request(row["url"], headers={"User-Agent": "landmark-research-reproduction/1.0"})
                with urlopen(request, timeout=90) as response, temp.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024**2)
                if digest(temp, row["algorithm"]) != row["checksum"]:
                    raise ValueError(f"Checksum mismatch: {path.name}")
                temp.replace(path)
                print(f"Retrieved {path.name}")
            finally:
                temp.unlink(missing_ok=True)


if __name__ == "__main__":
    fetch()
