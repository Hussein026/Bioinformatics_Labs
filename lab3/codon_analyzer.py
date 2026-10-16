"""
Codon frequency analyzer.
Reads a FASTA file (one virus/bacteria genome), slides a window of 3 over the
sequence, counts every codon, and shows the top 10 with their amino acids.

Usage: python codon_analyzer.py genome.fasta
Requires: pip install matplotlib
"""
import sys
import gzip
from collections import Counter

# RNA codon table (DNA codon -> RNA via T->U, then translate)
CODON_TABLE = {
    "UUU": "Phenylalanine", "UUC": "Phenylalanine", "UUA": "Leucine", "UUG": "Leucine",
    "CUU": "Leucine", "CUC": "Leucine", "CUA": "Leucine", "CUG": "Leucine",
    "AUU": "Isoleucine", "AUC": "Isoleucine", "AUA": "Isoleucine", "AUG": "Methionine",
    "GUU": "Valine", "GUC": "Valine", "GUA": "Valine", "GUG": "Valine",
    "UCU": "Serine", "UCC": "Serine", "UCA": "Serine", "UCG": "Serine",
    "CCU": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "ACU": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "GCU": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "UAU": "Tyrosine", "UAC": "Tyrosine", "UAA": "STOP", "UAG": "STOP",
    "CAU": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "AAU": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "GAU": "Aspartic acid", "GAC": "Aspartic acid", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "UGU": "Cysteine", "UGC": "Cysteine", "UGA": "STOP", "UGG": "Tryptophan",
    "CGU": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "AGU": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GGU": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
}

# Common food sources for each amino acid
FOODS = {
    "Alanine": "meat, fish, eggs, dairy, beans",
    "Arginine": "nuts, seeds, meat, fish, legumes",
    "Asparagine": "asparagus, potatoes, eggs, meat, whey",
    "Aspartic acid": "meat, sugar cane, asparagus, legumes",
    "Cysteine": "eggs, poultry, garlic, onions, broccoli",
    "Glutamic acid": "meat, fish, eggs, dairy, soy",
    "Glutamine": "beef, eggs, dairy, tofu, cabbage",
    "Glycine": "bone broth, gelatin, meat, fish, legumes",
    "Histidine": "meat, fish, dairy, nuts, whole grains",
    "Isoleucine": "chicken, fish, eggs, lentils, cheese",
    "Leucine": "beef, chicken, fish, eggs, soy, dairy",
    "Lysine": "red meat, fish, eggs, dairy, legumes",
    "Methionine": "eggs, fish, meat, Brazil nuts, sesame seeds",
    "Phenylalanine": "meat, fish, eggs, dairy, soy",
    "Proline": "meat, dairy, eggs, cabbage, asparagus",
    "Serine": "eggs, soy, meat, dairy, nuts",
    "Threonine": "meat, dairy, eggs, lentils, sesame seeds",
    "Tryptophan": "turkey, eggs, cheese, nuts, seeds, tofu",
    "Tyrosine": "cheese, soy, chicken, fish, almonds",
    "Valine": "meat, dairy, mushrooms, peanuts, soy",
    "STOP": "- (stop signal, not an amino acid)",
}


def read_fasta(path):
    """Return the genome sequence (uppercase) from a FASTA file."""
    seq = []
    opener = gzip.open if path.endswith('.gz') else open  # supports .fna.gz from NCBI
    with opener(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">") or not line:
                continue
            seq.append(line.upper())
    return "".join(seq)


def count_codons(seq):
    """Slide a window of 3 positions (step 1) and count each codon."""
    counts = Counter()
    for i in range(len(seq) - 2):
        codon = seq[i:i + 3]
        if set(codon) <= set("ACGT"):  # skip N or other ambiguous bases
            counts[codon] += 1
    return counts


def show_chart(top10):
    """Bar chart of the top 10 codons (opens a window and saves a PNG)."""
    import matplotlib.pyplot as plt

    labels = [f"{c}\n{CODON_TABLE[c.replace('T', 'U')]}" for c, _ in top10]
    values = [n for _, n in top10]
    plt.figure(figsize=(12, 6))
    bars = plt.bar(labels, values, color="steelblue")
    plt.bar_label(bars)
    plt.title("Top 10 codons in the genome")
    plt.xlabel("Codon (amino acid)")
    plt.ylabel("Count")
    plt.xticks(fontsize=8)
    plt.tight_layout()
    plt.savefig("top10_codons.png", dpi=150)
    plt.show()


def main():
    # Works with the VS Code Run button too: defaults to genome.fasta
    path = sys.argv[1] if len(sys.argv) > 1 else "genome.fasta"
    seq = read_fasta(path)
    counts = count_codons(seq)
    total = sum(counts.values())

    print(f"Genome length: {len(seq)} bases")
    print(f"Total codons counted: {total}\n")
    print(f"{'Rank':<5}{'DNA':<6}{'RNA':<6}{'Count':<9}{'Freq %':<9}Amino acid")
    print("-" * 50)

    top10 = counts.most_common(10)
    for rank, (codon, n) in enumerate(top10, 1):
        rna = codon.replace("T", "U")
        aa = CODON_TABLE[rna]
        print(f"{rank:<5}{codon:<6}{rna:<6}{n:<9}{n / total * 100:<9.2f}{aa}")

    print("\nFood sources for the amino acids of the top 10 codons:")
    seen = set()
    for codon, _ in top10:
        aa = CODON_TABLE[codon.replace("T", "U")]
        if aa not in seen:
            seen.add(aa)
            print(f"  {aa}: {FOODS[aa]}")

    show_chart(top10)


if __name__ == "__main__":
    main()
