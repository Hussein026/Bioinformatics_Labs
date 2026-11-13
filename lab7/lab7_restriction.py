"""
Bioinformatics - Lab 7 (part 2): restriction enzyme digestion + gel simulation

Input : 1. recognition site + cut position of each enzyme (table below)
        2. a DNA sequence (FASTA from NCBI)
Output: 1. number of cleavages for each enzyme
        2. cleavage positions and fragment lengths
        3. simulation of an electrophoresis gel (one lane per enzyme + ladder)

Positions are 0-based. A cut at position p means the DNA is cut between
base p-1 and base p, so fragment lengths are differences of cut positions.

Run:  python lab7_restriction.py
      python lab7_restriction.py --fasta sequence.fasta --circular
"""

import argparse
import math
import os
import urllib.request

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

DEFAULT_FASTA = "sequence.fasta"
DEFAULT_ACCESSION = "L09137"      # pUC19, used only if the FASTA is missing

# name: (recognition site 5'->3', number of bases from the site start to the cut)
#   EcoRI   G^AATTC      BamHI  G^GATCC     HindIII A^AGCTT
#   TaqI    T^CGA        HaeIII GG^CC (blunt)
ENZYMES = {
    "EcoRI":   ("GAATTC", 1),
    "BamHI":   ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI":    ("TCGA", 1),
    "HaeIII":  ("GGCC", 2),
}


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def download_fasta(accession, path):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={accession}&rettype=fasta&retmode=text")
    print(f"Downloading {accession} from NCBI ...")
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read().decode()
    with open(path, "w") as f:
        f.write(data)


def read_fasta(path):
    header, parts = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    break
                header = line[1:]
            else:
                parts.append(line.upper())
    seq = "".join(parts)
    return header, "".join(c for c in seq if c in BASES)


def find_cuts(seq, site, offset):
    """Cut positions of one enzyme. Overlapping sites are found too.
    All 5 enzymes here are palindromic (site == its reverse complement),
    so searching one strand is enough; other enzymes get both strands."""
    cuts = []
    pos = seq.find(site)
    while pos != -1:
        cuts.append(pos + offset)
        pos = seq.find(site, pos + 1)
    rc = reverse_complement(site)
    if rc != site:                               # non-palindromic enzyme
        pos = seq.find(rc)
        while pos != -1:
            cuts.append(pos + len(site) - offset)
            pos = seq.find(rc, pos + 1)
    return sorted(cuts)


def fragments(cuts, seq_len, circular):
    """Fragment lengths from the cut positions."""
    if not cuts:
        return [seq_len]
    if circular:
        if len(cuts) == 1:
            return [seq_len]                     # one cut opens the circle
        sizes = [b - a for a, b in zip(cuts, cuts[1:])]
        sizes.append(seq_len - cuts[-1] + cuts[0])
        return sizes
    edges = [0] + cuts + [seq_len]
    return [b - a for a, b in zip(edges, edges[1:]) if b - a > 0]


def digest(seq, circular):
    results = {}
    for name, (site, offset) in ENZYMES.items():
        cuts = find_cuts(seq, site, offset)
        results[name] = {"site": site, "cuts": cuts,
                         "fragments": fragments(cuts, len(seq), circular)}
    return results


def print_results(results, seq_len, circular):
    print(f"\nDNA is treated as {'circular' if circular else 'linear'}, {seq_len} bp\n")
    for name, r in results.items():
        print(f"{name}  (site {r['site']})")
        print(f"  number of cleavages : {len(r['cuts'])}")
        print(f"  cleavage positions  : {r['cuts'] if r['cuts'] else '-'}")
        print(f"  fragment lengths    : {r['fragments']}")
        print()


def make_ladder(max_size):
    ladder = [100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000,
              5000, 6000, 8000, 10000, 15000, 20000]
    top = [s for s in ladder if s < max_size * 1.5]
    return top if len(top) >= 5 else ladder[:8]


def draw_gel(results, seq_len, title, out_png, show):
    """Gel simulation: smaller fragments travel farther (log scale)."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lanes = [("Ladder", make_ladder(seq_len))]
    lanes.append(("Uncut", [seq_len]))
    lanes += [(n, r["fragments"]) for n, r in results.items()]

    all_sizes = [s for _, sizes in lanes for s in sizes]
    lo, hi = math.log10(max(min(all_sizes), 10) * 0.7), math.log10(max(all_sizes) * 1.4)

    def dist(size):                              # 0 = well (top), 1 = bottom
        return (hi - math.log10(size)) / (hi - lo)

    fig, ax = plt.subplots(figsize=(1.0 + 1.1 * len(lanes), 7))
    fig.patch.set_facecolor("#1b1b1f")
    ax.set_facecolor("#1b1b1f")
    for i, (name, sizes) in enumerate(lanes):
        ax.add_patch(plt.Rectangle((i - 0.35, -0.04), 0.7, 0.03, color="#555"))  # well
        # bands closer than 2% of the gel height run together: share one label
        groups = []
        for s in sorted(sizes, reverse=True):
            if groups and abs(dist(s) - dist(groups[-1][0])) < 0.02:
                groups[-1].append(s)
            else:
                groups.append([s])
        for g in groups:
            y = dist(g[0])
            strong = 0.9 if name != "Ladder" else 0.55
            ax.add_patch(plt.Rectangle((i - 0.33, y - 0.008), 0.66, 0.016,
                                       color="#e8f0ff", alpha=strong))
            label = " + ".join(str(x) for x in g)
            if name == "Ladder":
                ax.text(i - 0.4, y, label, color="#aab", fontsize=7, ha="right", va="center")
            else:
                ax.text(i + 0.4, y, label, color="#9fd", fontsize=7, ha="left", va="center")
        ax.text(i, -0.09, name, color="white", ha="center", fontsize=9)
    ax.set_xlim(-1, len(lanes))
    ax.set_ylim(1.05, -0.14)
    ax.axis("off")
    ax.set_title(title, color="white", fontsize=10, pad=12)
    ax.text(len(lanes) - 0.5, 1.03, "bp", color="#aab", fontsize=8, ha="right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, facecolor=fig.get_facecolor())
    print(f"Gel image saved to {out_png}")
    if show:
        plt.show()


def main():
    p = argparse.ArgumentParser(description="Restriction digestion + gel simulation")
    p.add_argument("--fasta", default=DEFAULT_FASTA)
    p.add_argument("--accession", default=DEFAULT_ACCESSION)
    p.add_argument("--circular", action="store_true",
                   help="treat the DNA as a circle (plasmids / bacterial genomes)")
    p.add_argument("--no-show", action="store_true", help="save the gel without opening a window")
    args = p.parse_args()

    if not os.path.exists(args.fasta):
        download_fasta(args.accession, args.fasta)
    header, seq = read_fasta(args.fasta)
    print(f"Sequence: {header}\nLength  : {len(seq)} bp")

    results = digest(seq, args.circular)
    print_results(results, len(seq), args.circular)
    draw_gel(results, len(seq), f"Simulated gel - {header[:50]}",
             "gel_simulation.png", show=not args.no_show)


if __name__ == "__main__":
    main()
