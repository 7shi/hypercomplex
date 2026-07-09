from octonion import triples

# Pure algebraic simulation (no matrices):
# L_1..L_7 are anticommuting symbols with L_i^2 = -1.
# Words of length >= 4 are folded to their complement (length <= 3)
# via the volume relation L_1 L_2 L_3 L_4 L_5 L_6 L_7 = -1.

def reduce_word(word):
    """Sort symbols by anticommutation, cancel L_i L_i = -1. Returns (sign, tuple)."""
    sign, w = 1, list(word)
    changed = True
    while changed:
        changed = False
        for i in range(len(w) - 1):
            if w[i] == w[i+1]:
                del w[i:i+2]; sign = -sign; changed = True; break
            if w[i] > w[i+1]:
                w[i], w[i+1] = w[i+1], w[i]; sign = -sign; changed = True; break
    return sign, tuple(w)

FULL = (1,2,3,4,5,6,7)

def fold(sign, word):
    """Replace a word of length >= 4 with its complement via L_1...L_7 = -1."""
    if len(word) < 4:
        return sign, word
    comp = tuple(x for x in FULL if x not in word)
    sigma, full = reduce_word(word + comp)   # word * comp = sigma * L_1...L_7 = -sigma
    assert full == FULL
    m = len(comp)
    inv_sign = (-1) ** (m * (m + 1) // 2)    # comp^{-1} = inv_sign * comp
    return sign * -sigma * inv_sign, comp    # word = -sigma * comp^{-1}

def fmt(sign, word):
    s = "-" if sign < 0 else " "
    return s + ("I" if not word else "".join(f"L{x}" for x in word))

# the 8 terms of P (as written in the article, with L_7 kept as a symbol)
p_terms = [(1, ())] + [(-1, t) for t in triples]

def multiply(prefix, label):
    print(f"--- 8 * {label} P ---")
    results = {}
    for s, t in p_terms:
        sign, word = reduce_word(prefix + t)
        fs, fw = fold(sign * s, word)
        note = "" if (fs, fw) == (sign * s, word) else f"   (= {fmt(sign * s, word)})"
        results[fw] = fs
        print(f"  {label} * {fmt(s, t):11s} = {fmt(fs, fw):11s}{note}")
    print()
    return results

a = multiply((1,2), "L1L2")
b = multiply((3,),  "L3  ")

print("same 8 terms with same signs:", a == b)

# the e_3 row of the classification table in the article
e3_row = {(): 0}  # placeholder replaced below
e3_row = {(3,): 1, (1,2): 1, (5,6): -1, (4,7): 1,
          (2,4,5): -1, (1,4,6): 1, (1,5,7): 1, (2,6,7): 1}
print("matches the e_3 row of the table :", a == e3_row)
