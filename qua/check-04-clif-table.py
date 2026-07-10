"""Checks for 04-clif-str-repr.md.

Verifies the complete classification table of Cl_{p,q}(R) (0 <= p,q <= 8)
in two independent ways:

1. Against the standard mod-8 classification (Bott periodicity): the type
   (R / 2R / C / H / 2H) depends only on p-q mod 8, and the matrix size is
   fixed by the total real dimension 2^(p+q).
2. By rebuilding the table with the article's own method: the series
   Cl_{p,q} (x) H = Cl_{q,p+2} and Cl_{p,q} (x) H' = Cl_{q+2,p}, starting
   from R, C, C', H, H', then filling the rest with
   Cl_{p+1,q+1} = Cl_{p,q} (x) H' = Cl_{p,q} (x) R(2),
   using the tensor rules F(m) (x) R(2) = F(2m), R (x) H = H,
   C (x) H = C(2), H (x) H = R(4).

Also checks the isomorphisms quoted in the tensor-rule section and the
8-periodicity section: H' (x) H = H(2), C (x) H = C (x) H' = C(2),
H (x) H = H' (x) H' = R(4), H(2) (x) H(2) = R(16), and the type-by-
(p-q mod 8) row stated at the end of the article.
"""

# an algebra "k F(n)" is represented as (k, F, n); k in {1, 2}
DIM = {"R": 1, "C": 2, "H": 4}

def parse(s):
    k = 1
    if s.startswith("2"):
        k, s = 2, s[1:]
    f = s[0]
    n = int(s[2:-1]) if "(" in s else 1
    return (k, f, n)

def show(a):
    k, f, n = a
    s = f if n == 1 else f"{f}({n})"
    return ("2" + s) if k == 2 else s

# the article's complete table (section 4.3), rows p = 0..8, columns q = 0..8
article = [row.split() for row in [
    "R     C     H     2H    H(2)  C(4)  R(8)  2R(8)  R(16)",
    "2R    R(2)  C(2)  H(2)  2H(2) H(4)  C(8)  R(16)  2R(16)",
    "R(2)  2R(2) R(4)  C(4)  H(4)  2H(4) H(8)  C(16)  R(32)",
    "C(2)  R(4)  2R(4) R(8)  C(8)  H(8)  2H(8) H(16)  C(32)",
    "H(2)  C(4)  R(8)  2R(8) R(16) C(16) H(16) 2H(16) H(32)",
    "2H(2) H(4)  C(8)  R(16) 2R(16) R(32) C(32) H(32)  2H(32)",
    "H(4)  2H(4) H(8)  C(16) R(32) 2R(32) R(64) C(64)  H(64)",
    "C(8)  H(8)  2H(8) H(16) C(32) R(64) 2R(64) R(128) C(128)",
    "R(16) C(16) H(16) 2H(16) H(32) C(64) R(128) 2R(128) R(256)",
]]

# --- check 1: standard mod-8 classification --------------------------------
# type of Cl_{p,q} by p-q mod 8 (convention: p generators square to +1)
MOD8 = {0: (1, "R"), 1: (2, "R"), 2: (1, "R"), 3: (1, "C"),
        4: (1, "H"), 5: (2, "H"), 6: (1, "H"), 7: (1, "C")}

def classify(p, q):
    k, f = MOD8[(p - q) % 8]
    n2 = 2 ** (p + q) // (k * DIM[f])  # n^2 from total real dimension
    n = round(n2 ** 0.5)
    assert n * n == n2
    return (k, f, n)

ok = all(parse(article[p][q]) == classify(p, q)
         for p in range(9) for q in range(9))
print("table matches the mod-8 classification:", ok)

# --- check 2: rebuild via the article's tensor-product method --------------
def tensor_r2(a):  # F(m) (x) R(2) = F(2m)
    k, f, n = a
    return (k, f, 2 * n)

def tensor_h(a):  # (x) H, using R(x)H=H, C(x)H=C(2), H(x)H=R(4)
    k, f, n = a
    return {"R": (k, "H", n), "C": (k, "C", 2 * n), "H": (k, "R", 4 * n)}[f]

table = {}
table[0, 0] = parse("R")     # Cl_{0,0}
table[0, 1] = parse("C")     # Cl_{0,1}
table[1, 0] = parse("2R")    # Cl_{1,0} = C' = R + R
table[0, 2] = parse("H")     # Cl_{0,2}
table[2, 0] = parse("R(2)")  # Cl_{2,0} = H' = R(2)

# extend the first row/column with (x)H : (p,q) -> (q,p+2)
# and (x)H' = (x)R(2) : (p,q) -> (q+2,p), as in section 3
for _ in range(4):
    for (p, q), a in list(table.items()):
        table[q, p + 2] = tensor_h(a)
        table[q + 2, p] = tensor_r2(a)

# fill the rest with Cl_{p+1,q+1} = Cl_{p,q} (x) H' = Cl_{p,q} (x) R(2)
for p in range(1, 9):
    for q in range(1, 9):
        if (p, q) not in table:
            table[p, q] = tensor_r2(table[p - 1, q - 1])

ok = all(table[p, q] == parse(article[p][q])
         for p in range(9) for q in range(9))
print("table matches the article's own construction:", ok)

# --- isomorphisms quoted in section 2 --------------------------------------
# H' (x) H = R(2) (x) H = H(2);  C (x) H = C (x) H' = C(2);  H (x) H = R(4)
print("H'(x)H = H(2):", tensor_h(parse("R(2)")) == parse("H(2)"))
print("C(x)H = C(x)H' = C(2):",
      tensor_h(parse("C")) == tensor_r2(parse("C")) == parse("C(2)"))
print("H(x)H = H'(x)H' = R(4):",
      tensor_h(parse("H")) == tensor_r2(tensor_r2(parse("R"))) == parse("R(4)"))

# --- 8-periodicity section --------------------------------------------------
# H(2) (x) H(2) = (H (x) H)(4) = R(16), so Cl_{p+8,q} = Cl_{p,q} (x) R(16)
k, f, n = tensor_h(parse("H(2)"))
print("H(2)(x)H(2) = R(16):", (k, f, 2 * n) == parse("R(16)"))

# type by p-q mod 8 as stated in the article
stated = "R 2R R C H 2H H C".split()
ok = all(classify(p, q)[:2] == parse(stated[(p - q) % 8])[:2]
         for p in range(9) for q in range(9))
print("type depends only on p-q mod 8 as stated:", ok)
