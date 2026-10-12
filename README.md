# fmindex_search

A compressed full text search index in pure Python, built from scratch with zero dependencies: linear time SA IS suffix array construction, the Burrows Wheeler transform, a wavelet matrix for rank queries, and FM index backward search with count, locate and extract.

## Why it matters

An FM index lets you search for any substring in a text in time proportional to the pattern length, not the text length, while storing the text in a form close to its compressed size. The same idea powers read aligners such as Bowtie and BWA in genomics. This project implements the whole stack so each layer can be read, tested and verified on its own.

## Layers

1. `fmindex/sais.py`: SA IS (Nong, Zhang, Chan) induced sorting suffix array in O(n).
2. `fmindex/wavelet.py`: a bit vector with block prefix rank, and a wavelet matrix giving O(log sigma) `access` and `rank` over any integer sequence.
3. `fmindex/index.py`: `FMIndex` with the BWT, the C array, LF mapping, backward search, sampled suffix array for `locate`, and `extract` to recover text from the index alone.

## Usage

```python
from fmindex import FMIndex

fm = FMIndex("mississippi", sa_sample=4)
fm.count("ssi")        # 2
fm.locate("ssi")       # [2, 5]
fm.extract(4, 5)       # b"issip"
fm.bwt()               # b"ipssm$pissii"
```

Command line:

```
python -m fmindex count  file.txt needle
python -m fmindex locate file.txt needle
python -m fmindex bwt    file.txt
```

`sa_sample` trades space for locate speed: only suffix array entries at text positions divisible by the rate are stored, and the rest are recovered by walking LF steps.

## Running the tests

```
python -m unittest discover -s tests
```

## Design and testing notes

- Text is treated as bytes, so UTF 8, binary data and arbitrary alphabets work. Symbols are shifted by one so that 0 is a unique smallest sentinel.
- The suffix array is checked against a naive sort of all suffixes on hundreds of random strings over alphabets of size 1, 2, 3 and 20, plus periodic strings which are the classic edge cases for induced sorting.
- The wavelet matrix is checked against list slicing for access and rank, including across 64 bit word boundaries.
- `FMIndex` count and locate are checked against a naive scan for several sampling rates. Extraction, binary input, non ASCII text, and absent patterns each have a test, and the BWT of banana is pinned to its known value.

## Limitations

This is an educational implementation. Rank uses Python integers and `bin().count`, so it is far slower than a C implementation, and the whole index is held in memory.

## License

MIT
