"""CLI: python -m fmindex {count,locate,bwt} FILE PATTERN"""
import argparse
import sys
from .index import FMIndex


def main(argv=None):
    p = argparse.ArgumentParser(prog="fmindex", description=__doc__)
    p.add_argument("cmd", choices=["count", "locate", "bwt"])
    p.add_argument("file")
    p.add_argument("pattern", nargs="?")
    p.add_argument("--sample", type=int, default=16, help="SA sampling rate")
    a = p.parse_args(argv)
    with open(a.file, "rb") as f:
        idx = FMIndex(f.read(), a.sample)
    if a.cmd == "bwt":
        sys.stdout.write(idx.bwt().decode("latin-1") + "\n")
        return 0
    if a.pattern is None:
        p.error("pattern required")
    if a.cmd == "count":
        print(idx.count(a.pattern))
    else:
        print("\n".join(map(str, idx.locate(a.pattern))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
