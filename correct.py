# Cell-barcode correction and UMI resolution with Hamming-distance-1 neighbour search.
from collections import Counter, defaultdict

BASES = "ACGT"


class Hamming1Index:
    """Find stored k-mers at Hamming distance <= 1 without scanning them all. Pigeonhole: two strings
    that differ in at most one position agree exactly on at least one half, so bucket by each half
    and only compare within the matching buckets. (Arcane builds on Fourway, a four-part refinement.)"""

    def __init__(self, kmers):
        self.kmers = list(kmers)
        self.half = len(self.kmers[0]) // 2
        self.left, self.right = defaultdict(list), defaultdict(list)
        for i, s in enumerate(self.kmers):
            self.left[s[:self.half]].append(i)
            self.right[s[self.half:]].append(i)

    def neighbours(self, q):
        cand = set(self.left.get(q[:self.half], ())) | set(self.right.get(q[self.half:], ()))
        return [self.kmers[i] for i in sorted(cand) if hamming(self.kmers[i], q) == 1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def correct_barcodes(observed, whitelist):
    """Exact whitelist hits stay. Others are corrected to a whitelist barcode at distance 1; if several
    qualify, the one with the most exact-match reads wins, and true ties are discarded."""
    wl = set(whitelist)
    exact_counts = Counter(b for b in observed if b in wl)
    index = Hamming1Index(sorted(exact_counts))  # only barcodes actually seen exactly
    fixed, status = [], []
    for b in observed:
        if b in wl:
            fixed.append(b); status.append("exact")
            continue
        cands = index.neighbours(b)
        if not cands:
            fixed.append(None); status.append("no whitelist neighbour")
            continue
        best = sorted(cands, key=lambda c: -exact_counts[c])
        if len(best) > 1 and exact_counts[best[0]] == exact_counts[best[1]]:
            fixed.append(None); status.append("ambiguous")
        else:
            fixed.append(best[0]); status.append("corrected")
    return fixed, status


def count_umis(umi_counts, method="directional"):
    """Molecules behind one (cell, gene) given {umi: reads}.
    unique      every distinct UMI is a molecule
    cluster     UMIs linked by Hamming-1 chains collapse to one
    directional UMI-tools: a absorbs b only if HD(a,b)=1 and reads(a) >= 2*reads(b) - 1"""
    umis = sorted(umi_counts, key=lambda u: -umi_counts[u])
    if method == "unique" or len(umis) < 2:
        return len(umis)
    index = Hamming1Index(umis)
    adj = {u: index.neighbours(u) for u in umis}
    seen, molecules = set(), 0
    for u in umis:  # highest count first
        if u in seen:
            continue
        molecules += 1
        stack = [u]
        seen.add(u)
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if b in seen:
                    continue
                if method == "cluster" or umi_counts[a] >= 2 * umi_counts[b] - 1:
                    seen.add(b)
                    stack.append(b)
    return molecules
