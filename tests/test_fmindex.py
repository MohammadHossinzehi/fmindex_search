import os, random, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from fmindex import suffix_array, WaveletMatrix, FMIndex


def naive_sa(b):
    return sorted(range(len(b)), key=lambda i: b[i:])


def naive_find(t, p):
    return [i for i in range(len(t) - len(p) + 1) if t[i:i + len(p)] == p]


class TestSA(unittest.TestCase):
    def test_classic(self):
        self.assertEqual(suffix_array("banana"), [5, 3, 1, 0, 4, 2])

    def test_random(self):
        rng = random.Random(1)
        for _ in range(400):
            n = rng.randint(1, 60)
            k = rng.choice([1, 2, 3, 20])
            b = bytes(rng.randrange(k) + 97 for _ in range(n))
            self.assertEqual(suffix_array(b), naive_sa(b), b)

    def test_periodic(self):
        for b in [b"a" * 50, b"ab" * 40, b"abc" * 30, b"aab" * 25]:
            self.assertEqual(suffix_array(b), naive_sa(b))


class TestWavelet(unittest.TestCase):
    def test_rank_access(self):
        rng = random.Random(2)
        seq = [rng.randrange(37) for _ in range(500)]
        wm = WaveletMatrix(seq, 37)
        for i in range(0, 500, 7):
            self.assertEqual(wm.access(i), seq[i])
        for c in range(37):
            for i in (0, 1, 63, 64, 65, 250, 500):
                self.assertEqual(wm.rank(c, i), seq[:i].count(c))


class TestFM(unittest.TestCase):
    def test_against_naive(self):
        rng = random.Random(3)
        for sample in (1, 4, 16):
            for _ in range(60):
                t = bytes(rng.randrange(3) + 97 for _ in range(rng.randint(0, 120)))
                fm = FMIndex(t, sample)
                for _ in range(15):
                    p = bytes(rng.randrange(3) + 97 for _ in range(rng.randint(1, 5)))
                    want = naive_find(t, p)
                    self.assertEqual(fm.count(p), len(want))
                    self.assertEqual(fm.locate(p), want)

    def test_banana_bwt(self):
        self.assertEqual(FMIndex("banana").bwt(), b"annb$aa")

    def test_extract(self):
        t = b"the quick brown fox jumps over the lazy dog"
        fm = FMIndex(t)
        for s in range(0, len(t), 5):
            for l in (0, 1, 7, 100):
                self.assertEqual(fm.extract(s, l), t[s:s + l])

    def test_binary_and_unicode(self):
        t = bytes(range(256)) * 3
        fm = FMIndex(t, 8)
        self.assertEqual(fm.locate(bytes([255, 0])), [255, 511])
        self.assertEqual(FMIndex("héllo wörld héllo").count("héllo"), 2)

    def test_absent(self):
        self.assertEqual(FMIndex("abc").count("zz"), 0)


if __name__ == "__main__":
    unittest.main()
