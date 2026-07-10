from octonion import L, P
import numpy as np

# 8P = I - L1L2L3 - L1L4L5 - L1L7L6 - L2L4L6 - L2L5L7 - L3L4L7 - L3L6L5
terms_8P = [
    (1, []),
    (-1, [1,2,3]),
    (-1, [1,4,5]),
    (-1, [1,7,6]),
    (-1, [2,4,6]),
    (-1, [2,5,7]),
    (-1, [3,4,7]),
    (-1, [3,6,5]),
]

def fmt(sign, idx):
    s = "+" if sign > 0 else "-"
    body = "I" if not idx else "".join(f"L_{k}" for k in idx)
    return f"{s}{body}"

def prod(idx):
    M = np.eye(8, dtype=int)
    for k in idx:
        M = M @ L[k]
    return M

# Repeated factors not at the front (e.g. -L5L1L4L5 in 8L_5P) are reduced by
# anticommuting the second occurrence next to the first, then applying
# L_i^2 = -I. Example:
#   L5L1L4L5 = -L5L1L5L4 = L5L5L1L4 = -L1L4   (two anticommutations cancel,
#                                                then L5L5=-I flips the sign)
#   so -L5L1L4L5 = -(-L1L4) = L1L4
# This only reduces terms with a repeated factor; a term whose remaining
# 3 letters are all distinct (e.g. -L2L1L4L5 in 8L_2P) is untouched.
def cancel_all(sign, word):
    # repeatedly bring a repeated index adjacent (via anticommutation) and
    # remove it with L_i^2 = -I: factor (-1)**(i2-i1) per pair
    word = list(word)
    while True:
        dup = next((v for v in set(word) if word.count(v) == 2), None)
        if dup is None:
            return sign, word
        i1 = word.index(dup)
        i2 = len(word) - 1 - word[::-1].index(dup)
        sign *= (-1) ** (i2 - i1)
        word = word[:i1] + word[i1+1:i2] + word[i2+1:]

# A leftover 4-distinct-letter word is folded to <=3 letters using
# L1L2L3L4L5L6 = L7:
#  - if the word contains L7, substitute L7 -> L1L2L3L4L5L6 in place, e.g.
#      -L1L2L5L7 = -L1L2L5(L1L2L3L4L5L6)
#    then cancel_all() pairs up L1,L2,L5 with their substituted twins,
#    leaving -L3L4L6.
#  - if the word has no L7, right-multiply by L7L6L5L4L3L2L1 = I (harmless),
#    e.g.
#      -L1L2L4L6 = -L1L2L4L6(L7L6L5L4L3L2L1)
#    then cancel_all() pairs up L1,L2,L4,L6 with their twins in the
#    appended block, leaving +L7L5L3.
def fold(sign, word):
    # reduce a 4-letter word (all distinct) to <=3 letters via L1L2L3L4L5L6=L7
    word = list(word)
    if 7 in word:
        i = word.index(7)
        word = word[:i] + [1,2,3,4,5,6] + word[i+1:]
    else:
        word = word + [7,6,5,4,3,2,1]  # = I
    return cancel_all(sign, word)

print("=== raw 8L_iP (left-multiply 8P by L_i, no reordering) + same-factor + fold reduction ===")
for i in range(1, 8):
    raw = []
    for sign, idx in terms_8P:
        if idx and idx[0] == i:
            raw.append((-sign, idx[1:]))          # L_i^2 = -I, adjacent cancel
        else:
            raw.append((sign, [i] + idx))

    reduced = []
    for sign, word in raw:
        sign, word = cancel_all(sign, word)       # same-factor reduction
        if len(word) == 4:
            sign, word = fold(sign, word)          # L7 = L1...L6 fold
        reduced.append((sign, word))

    line = " ".join(fmt(s, w) for s, w in reduced)
    print(f"8L_{i}P = {line}")

    total = sum(s * prod(w) for s, w in reduced)
    ok = np.array_equal(total, 8 * (L[i] @ P))
    print("  check:", ok)
