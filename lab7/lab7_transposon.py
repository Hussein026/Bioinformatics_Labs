"""
Bioinformatics - Lab 7: brute-force transposon finder

Idea (see the IS1 slide): a transposon has an inverted repeat (IR) at each end.
The second IR is the reverse complement of the first one.

Algorithm
  1. Read a DNA sequence from a FASTA file (downloaded from NCBI if missing).
  2. Brute force: generate EVERY combination of A, C, G, T of length k.
  3. For each combination, find all its positions in the sequence.
  4. For each position, look for the reverse complement of the combination
     further along the sequence. If found -> the two copies delimit a transposon.
  5. Report start index and length, then classify how transposons relate:
       - separate  : they do not touch each other
       - overlap   : they intersect only partially
       - nested    : one is completely inside another
  6. Draw a chart.

Indexes are 0-based (Python style). "end" is exclusive, so length = end - start.

Run:  python lab7_transposon.py
      python lab7_transposon.py --fasta my.fasta --min-len 10 --max-len 12
"""

import argparse
import itertools
import os
import time
import urllib.request
from dataclasses import dataclass

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

# ---------------------------------------------------------------- defaults ---
DEFAULT_FASTA = "sequence.fasta"
DEFAULT_ACCESSION = "L09137"      # pUC19 plasmid, used only if the FASTA is missing
DEFAULT_MIN_LEN = 10              # length of the inverted repeat (IR)
DEFAULT_MAX_LEN = 10              
DEFAULT_MIN_SPAN = 50             # smallest whole transposon (IR1 start -> IR2 end)
DEFAULT_MAX_SPAN = 1000           # largest whole transposon


@dataclass
class Transposon:
    ir: str            # the first inverted repeat (5' end)
    ir_rc: str         # its reverse complement (3' end)
    ir_len: int
    start: int         # index where IR1 begins
    ir2_start: int     # index where IR2 begins
    end: int           # index just after IR2 (exclusive)

    @property
    def length(self):
        return self.end - self.start


# ------------------------------------------------------------ sequence I/O ---
def reverse_complement(seq):
    """Reverse complement of a DNA string: complement each base, then reverse."""
    return seq.translate(COMPLEMENT)[::-1]


def download_fasta(accession, path):
    """Download a nucleotide FASTA record from NCBI (E-utilities)."""
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={accession}&rettype=fasta&retmode=text")
    print(f"Downloading {accession} from NCBI ...")
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read().decode()
    with open(path, "w") as f:
        f.write(data)


def read_fasta(path):
    """Return (header, sequence) of the first record; keeps only A/C/G/T."""
    header, parts = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:           # second record -> stop
                    break
                header = line[1:]
            else:
                parts.append(line.upper())
    seq = "".join(parts)
    return header, "".join(c for c in seq if c in BASES)


# ------------------------------------------------------------ brute force ---
def generate_combinations(k):
    """Every string of length k over A, C, G, T (4^k of them)."""
    for combo in itertools.product(BASES, repeat=k):
        yield "".join(combo)


def find_all(seq, pattern):
    """All start positions of pattern in seq (overlapping matches included)."""
    positions = []
    pos = seq.find(pattern)
    while pos != -1:
        positions.append(pos)
        pos = seq.find(pattern, pos + 1)
    return positions


def find_transposons(seq, min_len, max_len, min_span, max_span):
    hits = []
    for k in range(min_len, max_len + 1):
        for ir in generate_combinations(k):
            if ir not in seq:                      # quick reject
                continue
            first = find_all(seq, ir)
            rc = reverse_complement(ir)
            second = first if rc == ir else find_all(seq, rc)
            if not second:
                continue
            for i in first:
                for j in second:
                    if j < i + k:                  # IR2 must come after IR1
                        continue
                    span = j + k - i
                    if min_span <= span <= max_span:
                        hits.append(Transposon(ir, rc, k, i, j, j + k))
    return hits


def merge_stems(hits, seq):
    """One real inverted repeat of length L is found L-k+1 times when we search
    with a shorter k (shifted by one base each time). Merge those shifted hits
    into a single transposon with its full (maximal) inverted repeat."""
    by_diag = {}
    for h in hits:
        by_diag.setdefault((h.ir_len, h.start + h.ir2_start), []).append(h)
    merged = []
    for (k, _), group in by_diag.items():
        group.sort(key=lambda h: h.start)
        run = [group[0]]
        for h in group[1:] + [None]:
            if h is not None and h.start == run[-1].start + 1:
                run.append(h)
                continue
            first = run[0]                       # outermost hit of the run
            ir_len = k + len(run) - 1
            end = first.end
            merged.append(Transposon(seq[first.start:first.start + ir_len],
                                     seq[end - ir_len:end],
                                     ir_len, first.start, end - ir_len, end))
            run = [h]
    return merged


# -------------------------------------------------------------- scenarios ---
def relation(a, b):
    """Relation of two transposons (a starts first)."""
    if b.start >= a.end:
        return "separate"
    if b.end <= a.end:
        return "nested"       # b is completely inside a
    return "overlap"          # partial intersection


def classify(hits):
    """Return pair list [(a, b, relation)] and a per-transposon category."""
    hits.sort(key=lambda h: (h.start, -h.end))
    pairs, category = [], {id(h): "separate" for h in hits}
    for x in range(len(hits)):
        for y in range(x + 1, len(hits)):
            rel = relation(hits[x], hits[y])
            pairs.append((hits[x], hits[y], rel))
            if rel == "nested":
                category[id(hits[x])] = category[id(hits[y])] = "nested"
            elif rel == "overlap":
                for t in (hits[x], hits[y]):
                    if category[id(t)] != "nested":
                        category[id(t)] = "overlap"
    return pairs, category


# ------------------------------------------------------------------ output ---
def print_results(hits, pairs):
    print(f"\nTransposons found: {len(hits)}")
    ir5, ir3 = "IR (5')", "IR (3')"
    print(f"{'#':>3} {ir5:<14} {ir3:<14} {'IR len':>6} "
          f"{'start':>7} {'end':>7} {'length':>7}")
    for n, h in enumerate(hits, 1):
        print(f"{n:>3} {h.ir:<14} {h.ir_rc:<14} {h.ir_len:>6} "
              f"{h.start:>7} {h.end:>7} {h.length:>7}")

    names = {id(h): n for n, h in enumerate(hits, 1)}
    counts = {"separate": 0, "overlap": 0, "nested": 0}
    shown = {"overlap": [], "nested": []}
    for a, b, rel in pairs:
        counts[rel] += 1
        if rel != "separate":
            shown[rel].append(f"#{names[id(a)]} & #{names[id(b)]}")
    print("\nRelations between transposon pairs:")
    print(f"  completely separate : {counts['separate']}")
    print(f"  partial overlap     : {counts['overlap']}  {', '.join(shown['overlap'][:10])}")
    print(f"  one inside another  : {counts['nested']}  {', '.join(shown['nested'][:10])}")


def draw_chart(hits, category, seq_len, title, out_png, show):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    colors = {"separate": "#2a9d8f", "overlap": "#e9a23b", "nested": "#c8553d"}
    hits = hits[:40]
    fig, ax = plt.subplots(figsize=(11, 1.2 + 0.45 * max(len(hits), 3)))
    for row, h in enumerate(hits):
        c = colors[category[id(h)]]
        ax.barh(row, h.length, left=h.start, height=0.5, color=c, alpha=0.45)
        ax.barh(row, h.ir_len, left=h.start, height=0.5, color=c)          # IR 5'
        ax.barh(row, h.ir_len, left=h.ir2_start, height=0.5, color=c)      # IR 3'
        ax.text(h.end + seq_len * 0.005, row,
                f"#{row + 1}  start={h.start}  len={h.length}", va="center", fontsize=8)
    ax.set_xlim(0, seq_len * 1.18)
    ax.set_yticks(range(len(hits)))
    ax.set_yticklabels([f"#{i + 1}" for i in range(len(hits))])
    ax.invert_yaxis()
    ax.set_xlabel("position in sequence (bp)")
    ax.set_title(title)
    ax.legend(handles=[Patch(color=colors[k], label=k) for k in colors]
              + [Patch(color="gray", label="dark ends = inverted repeats")],
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"\nChart saved to {out_png}")
    if show:
        plt.show()


# -------------------------------------------------------------------- main ---
def main():
    p = argparse.ArgumentParser(description="Brute-force transposon finder")
    p.add_argument("--fasta", default=DEFAULT_FASTA)
    p.add_argument("--accession", default=DEFAULT_ACCESSION)
    p.add_argument("--min-len", type=int, default=DEFAULT_MIN_LEN)
    p.add_argument("--max-len", type=int, default=DEFAULT_MAX_LEN)
    p.add_argument("--min-span", type=int, default=DEFAULT_MIN_SPAN)
    p.add_argument("--max-span", type=int, default=DEFAULT_MAX_SPAN)
    p.add_argument("--all", action="store_true",
                   help="keep every raw hit (do not merge shifted copies of the same repeat)")
    p.add_argument("--no-show", action="store_true", help="save the chart without opening a window")
    args = p.parse_args()

    if not os.path.exists(args.fasta):
        download_fasta(args.accession, args.fasta)
    header, seq = read_fasta(args.fasta)
    print(f"Sequence: {header}\nLength  : {len(seq)} bp")
    print(f"Searching IR lengths {args.min_len}-{args.max_len} "
          f"({sum(4 ** k for k in range(args.min_len, args.max_len + 1)):,} combinations) ...")

    t0 = time.time()
    hits = find_transposons(seq, args.min_len, args.max_len, args.min_span, args.max_span)
    if not args.all:
        hits = merge_stems(hits, seq)
    print(f"Done in {time.time() - t0:.1f} s")

    pairs, category = classify(hits)
    print_results(hits, pairs)
    if hits:
        draw_chart(hits, category, len(seq), f"Transposons in {header[:60]}",
                   "transposons_chart.png", show=not args.no_show)


if __name__ == "__main__":
    main()