import math
import matplotlib.pyplot as plt


def validate(sequence):
    sequence = sequence.upper().strip()
    if not sequence or any(base not in "ATGC" for base in sequence):
        raise ValueError("Invalid DNA sequence. Use only A, T, G, C.")
    return sequence


def tm_simple(sequence):
    """Tm = 4(G+C) + 2(A+T)"""
    sequence = validate(sequence)
    gc = sequence.count("G") + sequence.count("C")
    at = sequence.count("A") + sequence.count("T")
    return 4 * gc + 2 * at


def tm_salt(sequence, na_conc):
    """Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length"""
    sequence = validate(sequence)
    length = len(sequence)
    gc = sequence.count("G") + sequence.count("C")
    percent_gc = (gc / length) * 100
    return 81.5 + 16.6 * math.log10(na_conc) + 0.41 * percent_gc - 600 / length


def make_sequence(length, gc_count):
    """Sequence with a chosen number of G/C bases, the rest A/T."""
    gc_part = ("GC" * length)[:gc_count]
    at_part = ("AT" * length)[:length - gc_count]
    return gc_part + at_part


# ---------------- Data ----------------
S = "TAACGCTA"
NA = 0.05  # mol/L (50 mM)
LENGTH = 8

tm1 = tm_simple(S)
tm2 = tm_salt(S, NA)

sequences = [make_sequence(LENGTH, k) for k in range(LENGTH + 1)]
g_counts = [s.count("G") for s in sequences]
c_counts = [s.count("C") for s in sequences]
simple_tms = [tm_simple(s) for s in sequences]
salt_tms = [tm_salt(s, NA) for s in sequences]

print("Sequence:", S)
print("Tm (simple formula):", tm1, "°C")
print("Tm (salt-adjusted formula):", round(tm2, 2), "°C")

# ---------------- Charts ----------------
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
ax1, ax2, ax3, ax4 = axes.flatten()

# Chart 1: base frequency of S
bases = ["A", "T", "G", "C"]
counts = [S.count(b) for b in bases]
ax1.bar(bases, counts, color=["green", "red", "orange", "blue"])
ax1.set_title(f"1) Base Frequency in {S}")
ax1.set_xlabel("Base")
ax1.set_ylabel("Frequency")
ax1.set_yticks(range(0, max(counts) + 1))
for i, c in enumerate(counts):
    ax1.text(i, c + 0.05, str(c), ha="center")

# Chart 2: Tm of S with both formulas
names = ["Simple", "Salt-adjusted"]
values = [tm1, tm2]
ax2.bar(names, values, color=["purple", "teal"])
ax2.set_title(f"2) Tm of {S} (two formulas)")
ax2.set_ylabel("Tm (°C)")
for i, v in enumerate(values):
    ax2.text(i, v + 0.3, f"{v:.2f}", ha="center")

# Chart 3: frequency of G and C in each sequence
x = range(len(sequences))
width = 0.4
ax3.bar([i - width / 2 for i in x], g_counts, width, label="G", color="orange")
ax3.bar([i + width / 2 for i in x], c_counts, width, label="C", color="blue")
ax3.set_title("3) Frequency of G and C in each sequence")
ax3.set_ylabel("Frequency")
ax3.set_xticks(list(x))
ax3.set_xticklabels(sequences, rotation=45, fontsize=8)
ax3.legend()

# Chart 4: Tm of each sequence
ax4.bar([i - width / 2 for i in x], simple_tms, width, label="Simple formula", color="purple")
ax4.bar([i + width / 2 for i in x], salt_tms, width, label="Salt-adjusted formula", color="teal")
ax4.set_title("4) Tm of each sequence (different G+C)")
ax4.set_ylabel("Tm (°C)")
ax4.set_xticks(list(x))
ax4.set_xticklabels(sequences, rotation=45, fontsize=8)
ax4.axhline(0, color="black", linewidth=0.8)
ax4.legend()

plt.tight_layout()
plt.show()