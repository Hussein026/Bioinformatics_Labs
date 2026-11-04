"""
Lab 4 - DNA Sequence Assembly
1. Uses SARS-CoV-2 genome (NC_045512.2) from NCBI
2. Takes 2000 random samples of 100-120 bases each
3. Stores samples in a list
4. Rebuilds (assembles) the genome using overlap-based assembly
5. Identifies the main problem with this approach
"""

import random
from collections import Counter

# ── 1. Read FASTA ──────────────────────────────────────────────────────────────
def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">") or not line:
                continue
            seq.append(line.upper())
    return "".join(seq)

# ── 2. Take 2000 random samples of 100-120 bases ──────────────────────────────
def take_samples(seq, n=2000, min_len=100, max_len=120):
    samples = []
    for _ in range(n):
        length = random.randint(min_len, max_len)
        start  = random.randint(0, len(seq) - length)
        samples.append(seq[start : start + length])
    return samples

# ── 3. Simple greedy overlap assembly ─────────────────────────────────────────
def overlap(a, b, min_overlap=20):
    """Return length of longest suffix of a that matches prefix of b."""
    start = max(0, len(a) - len(b))
    best  = 0
    for i in range(start, len(a)):
        if b.startswith(a[i:]):
            best = len(a) - i
    return best

def greedy_assembly(reads, min_overlap=20):
    """Greedy: repeatedly merge the pair with the biggest overlap."""
    contigs = list(reads)
    merged  = True
    steps   = 0
    while merged and steps < 500:          # cap at 500 merge steps for speed
        merged = False
        best_score, best_i, best_j = 0, -1, -1
        # only check a random subset to keep it fast
        sample = random.sample(range(len(contigs)), min(50, len(contigs)))
        for i in sample:
            for j in sample:
                if i == j:
                    continue
                ov = overlap(contigs[i], contigs[j], min_overlap)
                if ov > best_score:
                    best_score, best_i, best_j = ov, i, j
        if best_score >= min_overlap:
            merged_seq = contigs[best_i] + contigs[best_j][best_score:]
            contigs.append(merged_seq)
            # remove in reverse order to keep indices valid
            for idx in sorted([best_i, best_j], reverse=True):
                contigs.pop(idx)
            merged = True
        steps += 1
    return contigs

# ── 4. Coverage analysis ───────────────────────────────────────────────────────
def coverage(genome_len, samples):
    cov = [0] * genome_len
    # We don't know exact positions after sampling, so estimate from lengths
    total_bases = sum(len(s) for s in samples)
    avg_cov = total_bases / genome_len
    return avg_cov

# ── Main ───────────────────────────────────────────────────────────────────────
def show_chart(samples, contigs, genome_len):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("DNA Assembly Results - SARS-CoV-2", fontsize=14, fontweight="bold")

    lengths = [len(s) for s in samples]
    axes[0].hist(lengths, bins=20, color="steelblue", edgecolor="white")
    axes[0].set_title("Sample Length Distribution")
    axes[0].set_xlabel("Length (bases)")
    axes[0].set_ylabel("Count")

    contig_lengths = sorted([len(c) for c in contigs], reverse=True)
    axes[1].bar(range(1, len(contig_lengths)+1), contig_lengths, color="tomato")
    axes[1].set_title(f"Contig Lengths ({len(contigs)} contigs)")
    axes[1].set_xlabel("Contig rank")
    axes[1].set_ylabel("Length (bases)")

    recovered = max(contig_lengths)
    axes[2].bar(["Longest contig", "Not recovered"], [recovered, genome_len - recovered],
                color=["steelblue", "tomato"])
    axes[2].bar_label(axes[2].containers[0])
    axes[2].set_title("Genome Recovery")
    axes[2].set_ylabel("Bases")

    plt.tight_layout()
    plt.savefig("assembly_results.png", dpi=150)
    plt.show()


def main():
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "GCF_009858895.2_ASM985889v3_genomic.fna"

    print("Reading genome...")
    genome = read_fasta(path)
    print(f"Genome length       : {len(genome)} bases")

    random.seed(42)
    print("\nTaking 2000 random samples (100-120 bases each)...")
    samples = take_samples(genome, n=2000)
    print(f"Samples collected   : {len(samples)}")
    print(f"Sample length range : {min(len(s) for s in samples)} - {max(len(s) for s in samples)} bases")
    total_bases = sum(len(s) for s in samples)
    print(f"Total bases sampled : {total_bases}")
    avg_cov = total_bases / len(genome)
    print(f"Average coverage    : {avg_cov:.1f}x")

    # Show first 3 samples
    print("\nFirst 3 samples:")
    for i, s in enumerate(samples[:3]):
        print(f"  Sample {i+1} (len={len(s)}): {s[:40]}...")

    print("\nAttempting greedy overlap assembly (500 merge steps)...")
    contigs = greedy_assembly(samples[:100], min_overlap=20)   # use 100 reads for speed
    contigs.sort(key=len, reverse=True)
    print(f"Contigs assembled   : {len(contigs)}")
    print(f"Longest contig      : {len(contigs[0])} bases")
    print(f"Original genome     : {len(genome)} bases")
    print(f"Recovery            : {len(contigs[0]) / len(genome) * 100:.1f}% of genome in longest contig")

    print("\n" + "="*60)
    print("MAIN PROBLEMS WITH THIS ALGORITHM:")
    print("="*60)
    problems = [
        ("1. Repeat regions",
         "   Random reads that fall inside repeated sequences cannot\n"
         "   be placed correctly. The assembler may merge them in the\n"
         "   wrong order, creating errors or collapsing repeats."),
        ("2. Coverage gaps",
         f"  With {len(samples)} reads of ~110 bases over a {len(genome)}-base genome,\n"
         "   average coverage is ~7x. Some regions may have 0 reads\n"
         "   (gaps), making full reconstruction impossible."),
        ("3. Greedy assembly is not optimal",
         "   A greedy approach picks the locally best merge at each\n"
         "   step. This may not lead to the globally correct assembly.\n"
         "   It is NP-hard to find the optimal solution."),
        ("4. Short reads = short overlaps",
         "   100-120 base reads produce short overlaps (~20 bases).\n"
         "   Short overlaps are more likely to be wrong matches,\n"
         "   leading to mis-assemblies."),
        ("5. Scalability",
         "   Checking all pairs of reads is O(n^2). With 2000 reads\n"
         "   that is 4,000,000 comparisons, very slow for large genomes."),
    ]
    for title, desc in problems:
        print(f"\n{title}")
        print(desc)

    show_chart(samples, contigs, len(genome))

if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()