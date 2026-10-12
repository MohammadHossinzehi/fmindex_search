"""Suffix array construction by induced sorting (SA-IS), O(n) time.

Reference: Nong, Zhang, Chan, "Two Efficient Algorithms for Linear Time
Suffix Array Construction" (2011). Input is a list of ints in [0, K).
The last symbol must be a unique, smallest sentinel (0).
"""
from __future__ import annotations


def _bucket_bounds(s, K, end):
    sizes = [0] * K
    for c in s:
        sizes[c] += 1
    out, total = [0] * K, 0
    for i in range(K):
        total += sizes[i]
        out[i] = total if end else total - sizes[i]
    return out


def _sais(s, K):
    n = len(s)
    if n == 1:
        return [0]
    if n == 2:
        return [1, 0]
    # classify: True = S-type
    t = [False] * n
    t[n - 1] = True
    for i in range(n - 2, -1, -1):
        t[i] = s[i] < s[i + 1] or (s[i] == s[i + 1] and t[i + 1])

    def is_lms(i):
        return i > 0 and t[i] and not t[i - 1]

    def induce(sa, lms_order):
        for i in range(n):
            sa[i] = -1
        ends = _bucket_bounds(s, K, True)
        for i in reversed(lms_order):
            ends[s[i]] -= 1
            sa[ends[s[i]]] = i
        starts = _bucket_bounds(s, K, False)
        for i in range(n):
            j = sa[i] - 1
            if sa[i] > 0 and not t[j]:
                sa[starts[s[j]]] = j
                starts[s[j]] += 1
        ends = _bucket_bounds(s, K, True)
        for i in range(n - 1, -1, -1):
            j = sa[i] - 1
            if sa[i] > 0 and t[j]:
                ends[s[j]] -= 1
                sa[ends[s[j]]] = j

    sa = [-1] * n
    lms = [i for i in range(1, n) if is_lms(i)]
    induce(sa, lms)

    # name LMS substrings
    sorted_lms = [i for i in sa if is_lms(i)]
    names = [-1] * n
    name, prev = 0, -1
    for i in sorted_lms:
        diff = prev == -1
        if not diff:
            for d in range(n):
                a, b = i + d, prev + d
                if s[a] != s[b] or t[a] != t[b]:
                    diff = True
                    break
                if d > 0 and (is_lms(a) or is_lms(b)):
                    break
        if diff:
            name += 1
        prev = i
        names[i] = name - 1
    reduced = [names[i] for i in lms]
    if name < len(lms):
        rsa = _sais(reduced, name)
    else:
        rsa = [0] * len(lms)
        for i, c in enumerate(reduced):
            rsa[c] = i
    induce(sa, [lms[i] for i in rsa])
    return sa


def suffix_array(data: bytes | str) -> list[int]:
    """Suffix array of ``data`` (sentinel handled internally, not returned)."""
    if isinstance(data, str):
        data = data.encode("utf-8")
    s = [b + 1 for b in data] + [0]
    return _sais(s, 257)[1:]
