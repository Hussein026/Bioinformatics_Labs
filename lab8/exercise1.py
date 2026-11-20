import matplotlib.pyplot as plt


def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
    n, m = len(s1), len(s2)

    # Initialize scoring matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        dp[i][0] = i * gap
    for j in range(m + 1):
        dp[0][j] = j * gap

    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if s1[i-1] == s2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + score,   # diagonal
                dp[i-1][j]   + gap,     # up (gap in s2)
                dp[i][j-1]   + gap      # left (gap in s1)
            )

    # Traceback
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = match if s1[i-1] == s2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + score:
                aligned1.append(s1[i-1])
                aligned2.append(s2[j-1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            aligned1.append(s1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(s2[j-1])
            j -= 1

    return (''.join(reversed(aligned1)),
            ''.join(reversed(aligned2)),
            dp)


def print_matrix(dp, s1, s2):
    header = " " * 8 + "  ".join(f"{c:>3}" for c in s2)
    print(header)
    for i, row in enumerate(dp):
        label = s1[i-1] if i > 0 else " "
        print(f"{label:2} " + "  ".join(f"{v:3}" for v in row))


def plot_matrix(dp, s1, s2, score):
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(dp, cmap="viridis")

    ax.set_xticks(range(len(s2) + 1))
    ax.set_xticklabels(["-"] + list(s2))
    ax.set_yticks(range(len(s1) + 1))
    ax.set_yticklabels(["-"] + list(s1))
    ax.xaxis.tick_top()

    for i in range(len(dp)):
        for j in range(len(dp[0])):
            ax.text(j, i, dp[i][j], ha="center", va="center",
                    color="white", fontsize=8)

    ax.set_title(f"Needleman-Wunsch Matrix (score = {score})", pad=30)
    fig.colorbar(im, label="Score")
    plt.tight_layout()
    plt.savefig("matrix.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    s1 = "ATCCGCCTAATAC"
    s2 = "CGAAATT"

    print(f"Sequence 1: {s1}")
    print(f"Sequence 2: {s2}")
    print("Scoring : match=+1, mismatch=-1, gap=-2\n")

    a1, a2, dp = needleman_wunsch(s1, s2)
    score = dp[len(s1)][len(s2)]

    print("Scoring Matrix:")
    print_matrix(dp, s1, s2)

    match_line = ''.join('|' if a == b else ' ' for a, b in zip(a1, a2))

    print("\nAlignment:")
    print(f"  {a1}")
    print(f"  {match_line}")
    print(f"  {a2}")
    print(f"\nAlignment score: {score}")

    plot_matrix(dp, s1, s2, score)