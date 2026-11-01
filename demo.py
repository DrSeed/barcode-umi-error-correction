# Self-contained demo: simulated droplet reads with barcode and UMI errors, resolved three ways.
import os
import time
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from correct import BASES, Hamming1Index, correct_barcodes, count_umis

rng = np.random.default_rng(5)


ASCII = np.frombuffer(BASES.encode(), dtype=np.uint8)


def to_str(codes):
    """Rows of base codes (0-3) -> Python strings, vectorised."""
    k = codes.shape[1]
    return ASCII[codes].view(f"S{k}").ravel().astype(f"U{k}").tolist()


def mutate(codes, rate):
    """Substitute each base with probability `rate`, always to a different base."""
    hit = rng.uniform(size=codes.shape) < rate
    return np.where(hit, (codes + rng.integers(1, 4, codes.shape)) % 4, codes)


def simulate(n_cells=100, n_genes=40, whitelist_size=20000, umi_len=10, bc_err=0.003, pcr_err=0.02):
    wl_codes = np.unique(rng.integers(0, 4, (whitelist_size, 16)), axis=0)
    whitelist = to_str(wl_codes)
    cell_codes = wl_codes[rng.choice(len(wl_codes), n_cells, replace=False)]
    gene_rate = rng.lognormal(1.5, 1.2, n_genes)
    gene_rate[0] = 1500  # one very highly expressed gene (think mitochondrial)
    size = rng.lognormal(0, 0.3, n_cells)
    n_mol = rng.poisson(np.outer(size, gene_rate))           # cells x genes, the truth
    mol_cell = np.repeat(np.arange(n_cells), n_mol.sum(1))
    mol_gene = np.concatenate([np.repeat(np.arange(n_genes), row) for row in n_mol])
    mol_umi = rng.integers(0, 4, (len(mol_cell), umi_len))
    dup = 1 + rng.poisson(2, len(mol_cell))                   # PCR duplicates per molecule
    r_mol = np.repeat(np.arange(len(mol_cell)), dup)
    bc = to_str(mutate(cell_codes[mol_cell[r_mol]], bc_err))  # sequencing errors in the barcode
    umi = to_str(mutate(mol_umi[r_mol], pcr_err / umi_len))   # PCR/sequencing errors in the UMI
    cells = to_str(cell_codes)
    truth = {(cells[c], g): int(n_mol[c, g]) for c in range(n_cells) for g in range(n_genes)}
    reads = list(zip(bc, mol_gene[r_mol].tolist(), umi))
    return whitelist, truth, reads


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    whitelist, truth, reads = simulate()
    true_cb = {cb for cb, _ in truth}

    fixed, status = correct_barcodes([r[0] for r in reads], whitelist)
    # speed: neighbours of 200 erroneous barcodes in the full whitelist, with and without the index
    wl = set(whitelist)
    sample = list(dict.fromkeys(r[0] for r in reads if r[0] not in wl))[:200]
    index = Hamming1Index(whitelist)
    t = time.time()
    hits_idx = [sorted(index.neighbours(b)) for b in sample]
    t_index = time.time() - t
    wl_arr = np.array([list(w) for w in whitelist])
    t = time.time()  # brute force, vectorised with NumPy: compare each query with every whitelist entry
    hits_brute = [sorted(np.array(whitelist)[(wl_arr != np.array(list(b))).sum(1) == 1]) for b in sample]
    t_brute = time.time() - t
    assert hits_idx == hits_brute

    groups = defaultdict(Counter)
    for (cb, g, umi), fx in zip(reads, fixed):
        if fx is not None:
            groups[(fx, g)][umi] += 1
    rows = []
    for key, n_true in truth.items():
        if n_true == 0:
            continue
        uc = groups.get(key, Counter())
        rows.append({"cell": key[0], "gene": key[1], "true": n_true,
                     **{m: count_umis(uc, m) for m in ["unique", "cluster", "directional"]}})
    est = pd.DataFrame(rows)
    est.to_csv("results/umi_counts.csv", index=False)

    st = Counter(status)
    miscorrected = sum(1 for (cb, _, _), fx, s in zip(reads, fixed, status) if s == "corrected" and fx not in true_cb)
    summary = [{"metric": f"reads: {k}", "value": round(100 * v / len(reads), 2), "unit": "% of reads"} for k, v in st.items()]
    summary.append({"metric": "corrected reads sent to a non-cell barcode", "value": miscorrected, "unit": "reads"})
    summary.append({"metric": "Hamming-1 lookup of 200 barcodes in the whitelist, half-split index", "value": round(t_index, 4), "unit": "seconds"})
    summary.append({"metric": "Hamming-1 lookup of 200 barcodes in the whitelist, brute force", "value": round(t_brute, 2), "unit": "seconds"})
    for m in ["unique", "cluster", "directional"]:
        rel = (est[m] - est["true"]) / est["true"]
        summary.append({"metric": f"UMI method {m}: median relative error", "value": round(100 * rel.median(), 1), "unit": "%"})
        hi = est["gene"] == 0
        summary.append({"metric": f"UMI method {m}: relative error, top gene", "value": round(100 * rel[hi].median(), 1), "unit": "%"})
    table = pd.DataFrame(summary)
    table.to_csv("results/summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    labels = ["corrected", "ambiguous", "no whitelist neighbour"]
    vals = [100 * st.get(k, 0) / len(reads) for k in labels]
    axes[0].bar(labels, vals, color=["#55a868", "#dd8452", "#c44e52"])
    for i, v in enumerate(vals):
        axes[0].text(i, v + 0.05, f"{v:.2f}%", ha="center")
    axes[0].set_ylabel("% of reads")
    axes[0].set_title(f"Cell barcodes ({100 * st['exact'] / len(reads):.1f}% exact, not shown): most errors are rescued\n"
                      f"Hamming-1 search in {len(whitelist):,} barcodes: {1000 * t_index:.0f} ms indexed vs {1000 * t_brute:.0f} ms NumPy brute force",
                      fontsize=10)

    colours = {"unique": "#c44e52", "cluster": "#dd8452", "directional": "#4c72b0"}
    bins = [0, 5, 20, 80, 300, 5000]
    est["bin"] = pd.cut(est["true"], bins)
    for k, (m, c) in enumerate(colours.items()):
        rel = 100 * (est[m] - est["true"]) / est["true"]
        med = rel.groupby(est["bin"], observed=True).median()
        axes[1].plot(range(len(med)), med.values, "o-", color=c, label=m)
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].set_xticks(range(len(bins) - 1), [f"{a + 1}-{b}" for a, b in zip(bins[:-1], bins[1:])])
    axes[1].set_xlabel("true molecules per cell and gene")
    axes[1].set_ylabel("median error in molecule count (%)")
    axes[1].set_title("Unique UMIs over-count; plain clustering\nunder-counts busy genes; directional holds up")
    axes[1].legend()

    hi = est[est["gene"] == 0]
    for m, c in colours.items():
        axes[2].scatter(hi["true"], hi[m], s=10, color=c, alpha=0.7, label=m)
    lim = [hi["true"].min() * 0.9, hi[["unique", "true"]].max().max() * 1.05]
    axes[2].plot(lim, lim, "k--", lw=0.8)
    axes[2].set_xlabel("true molecules (highly expressed gene)")
    axes[2].set_ylabel("estimated molecules")
    axes[2].set_title("Collisions: with 10-nt UMIs, a busy gene's real\nmolecules often sit 1 mismatch apart")
    axes[2].legend()
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
