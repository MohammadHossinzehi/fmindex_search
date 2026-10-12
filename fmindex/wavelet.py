"""Wavelet matrix with O(log sigma) rank/access over an integer sequence."""
from __future__ import annotations
from bisect import bisect_right


class BitVector:
    """Bit vector with O(1) rank via per-64-bit-block prefix counts."""

    def __init__(self, bits):
        self.n = len(bits)
        self.words = []
        self.cum = [0]
        for i in range(0, self.n, 64):
            w = 0
            for j, b in enumerate(bits[i:i + 64]):
                if b:
                    w |= 1 << j
            self.words.append(w)
            self.cum.append(self.cum[-1] + bin(w).count("1"))
        self.ones = self.cum[-1]

    def rank1(self, i):
        """Number of 1 bits in [0, i)."""
        q, r = divmod(i, 64)
        if q >= len(self.words):
            return self.ones
        return self.cum[q] + bin(self.words[q] & ((1 << r) - 1)).count("1")

    def rank0(self, i):
        return i - self.rank1(i)

    def get(self, i):
        return (self.words[i >> 6] >> (i & 63)) & 1


class WaveletMatrix:
    def __init__(self, seq, sigma=None):
        seq = list(seq)
        self.n = len(seq)
        sigma = sigma or (max(seq) + 1 if seq else 1)
        self.levels = max(1, (sigma - 1).bit_length())
        self.bv, self.zeros = [], []
        cur = seq
        for lvl in range(self.levels - 1, -1, -1):
            bits = [(x >> lvl) & 1 for x in cur]
            bv = BitVector(bits)
            self.bv.append(bv)
            self.zeros.append(self.n - bv.ones)
            cur = [x for x, b in zip(cur, bits) if not b] + \
                  [x for x, b in zip(cur, bits) if b]

    def access(self, i):
        v = 0
        for k, bv in enumerate(self.bv):
            b = bv.get(i)
            v = (v << 1) | b
            i = bv.rank1(i) + self.zeros[k] if b else bv.rank0(i)
        return v

    def rank(self, c, i):
        """Occurrences of symbol c in seq[0:i]."""
        lo, hi = 0, i
        for k, bv in enumerate(self.bv):
            b = (c >> (self.levels - 1 - k)) & 1
            if b:
                lo = bv.rank1(lo) + self.zeros[k]
                hi = bv.rank1(hi) + self.zeros[k]
            else:
                lo, hi = bv.rank0(lo), bv.rank0(hi)
        return hi - lo
