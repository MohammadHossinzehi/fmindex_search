"""FM-index: backward search over the Burrows-Wheeler transform."""
from __future__ import annotations
import pickle
from .sais import suffix_array
from .wavelet import WaveletMatrix


class FMIndex:
    def __init__(self, text: bytes | str, sa_sample: int = 16):
        if isinstance(text, str):
            text = text.encode("utf-8")
        if sa_sample < 1:
            raise ValueError("sa_sample must be >= 1")
        self.n = len(text)
        self.sample = sa_sample
        sa = [self.n] + suffix_array(text)          # include sentinel suffix
        # BWT symbols: byte+1, sentinel = 0
        bwt = [(text[p - 1] + 1) if p else 0 for p in sa]
        self.wm = WaveletMatrix(bwt, 257)
        counts = [0] * 258
        for c in bwt:
            counts[c + 1] += 1
        for i in range(1, 258):
            counts[i] += counts[i - 1]
        self.C = counts                              # C[c] = #symbols < c
        self.sampled = {i: p for i, p in enumerate(sa) if p % sa_sample == 0}
        self.sentinel_row = bwt.index(0)

    # LF mapping
    def _lf(self, i):
        c = self.wm.access(i)
        return self.C[c] + self.wm.rank(c, i)

    def _range(self, pattern: bytes):
        lo, hi = 0, self.n + 1
        for b in reversed(pattern):
            c = b + 1
            lo = self.C[c] + self.wm.rank(c, lo)
            hi = self.C[c] + self.wm.rank(c, hi)
            if lo >= hi:
                return 0, 0
        return lo, hi

    def count(self, pattern: bytes | str) -> int:
        if isinstance(pattern, str):
            pattern = pattern.encode("utf-8")
        if not pattern:
            return self.n + 1
        lo, hi = self._range(pattern)
        return hi - lo

    def _locate_row(self, i):
        steps = 0
        while i not in self.sampled:
            i = self._lf(i)
            steps += 1
        return self.sampled[i] + steps

    def locate(self, pattern: bytes | str) -> list[int]:
        if isinstance(pattern, str):
            pattern = pattern.encode("utf-8")
        if not pattern:
            return list(range(self.n + 1))
        lo, hi = self._range(pattern)
        return sorted(self._locate_row(i) for i in range(lo, hi))

    def extract(self, start: int, length: int) -> bytes:
        """Recover text[start:start+length] from the index alone."""
        end = min(self.n, start + length)
        if start < 0 or start > end:
            raise ValueError("bad range")
        # find row whose suffix starts at `end`: walk from sentinel row
        # (row 0 is the sentinel suffix at position n).
        i, pos, out = 0, self.n, []
        while pos > start:
            c = self.wm.access(i)
            pos -= 1
            if pos < end:
                out.append(c - 1)
            i = self.C[c] + self.wm.rank(c, i)
        return bytes(reversed(out))

    def bwt(self) -> bytes:
        return bytes(0x24 if (c := self.wm.access(i)) == 0 else c - 1
                     for i in range(self.n + 1))

    def save(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(path):
        with open(path, "rb") as f:
            return pickle.load(f)
