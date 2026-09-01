from sage.all import *
import time
import csv

# =====================================
# EKSPERYMENT
# =====================================

def experiment(bits=32):

    print("\n==== bits =", bits, "====")

    # =====================================
    # LOSOWE p
    # =====================================

    while True:
        p = random_prime(2**bits - 1, lbound=2**(bits - 1))
        if p % 4 == 3:
            break

    F = GF(p)
    E = EllipticCurve(F, [1, 0])
    N = E.cardinality()

    print("Group order N =", N)

    # =====================================
    # DUŻY CZYNNIK PIERWSZY
    # =====================================

    m = max([f[0] for f in factor(N)])

    print("Largest prime factor m =", m)
    print("bitlength(m) =", m.nbits())

    # =====================================
    # PUNKT RZĘDU m
    # =====================================

    P = (N // m) * E.random_point()

    while P.is_zero() or P.order() != m:
        P = (N // m) * E.random_point()

    # =====================================
    # TAJNY LOGARYTM
    # =====================================

    x = ZZ.random_element(1, m)
    R = x * P

    # =====================================
    # PRZED REDUKCJĄ MOV
    # =====================================

    start = time.time()
    x_before_mov = discrete_log(R, P, operation='+')
    time_before_mov = time.time() - start

    # =====================================
    # EMBEDDING DEGREE
    # =====================================

    k = Mod(p, m).multiplicative_order()
    print("embedding degree =", k)

    # =====================================
    # PRZEJŚCIE DO F_{p^k}
    # =====================================

    Ek = E.base_extend(GF(p**k))

    Pk = Ek(P)
    Rk = Ek(R)

    # =====================================
    # WYBÓR Q
    # =====================================

    while True:
        Q = Ek.random_point()
        alpha = Pk.tate_pairing(Q, m, k)
        if alpha.multiplicative_order() == m:
            break

    beta = Rk.tate_pairing(Q, m, k)

    # =====================================
    # PO REDUKCJI MOV
    # =====================================

    start = time.time()
    x_after_mov = discrete_log(beta, alpha)
    time_after_mov = time.time() - start

    print("time before MOV =", time_before_mov)
    print("time after MOV =", time_after_mov)

    return {
        "m": int(m),
        "m_bits": int(m.nbits()),
        "p": int(p),
        "p_bits": int(p.nbits()),
        "a": int(1),
        "b": int(0),
        "N": int(N),
        "embedding_degree": int(k),
        "time_before_mov": float(time_before_mov),
        "time_after_mov": float(time_after_mov)
    }


# =====================================
# ZBIERANIE DANYCH
# =====================================

all_data = []

for b in [28, 32, 36, 40]:

    for _ in range(50):

        try:

            res = experiment(bits=b)
            all_data.append(res)

        except Exception as e:

            print("ERROR:", e)

# =====================================
# FILTROWANIE
# =====================================

MIN_M = 10**5

filtered_data = [d for d in all_data if d["m"] > MIN_M]

print(f"\nPo filtracji m > {MIN_M}:")
print("liczba punktów =", len(filtered_data))

filtered_data.sort(key=lambda d: d["m"])

# =====================================
# DATASET
# =====================================

with open(
    "mov_dataset.csv",
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "m",
            "m_bits",
            "p",
            "p_bits",
            "a",
            "b",
            "N"
        ]
    )

    writer.writeheader()

    for row in filtered_data:

        writer.writerow({
            "m": row["m"],
            "m_bits": row["m_bits"],
            "p": row["p"],
            "p_bits": row["p_bits"],
            "a": row["a"],
            "b": row["b"],
            "N": row["N"]
        })

print(f"Zapisano {len(filtered_data)} rekordów do mov_dataset.csv")

# =====================================
# WYNIKI EKSPERYMENTÓW
# =====================================

with open(
    "mov_results.csv",
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "m",
            "m_bits",
            "p",
            "p_bits",
            "a",
            "b",
            "N",
            "embedding_degree",
            "time_before_mov",
            "time_after_mov"
        ]
    )

    writer.writeheader()

    for row in filtered_data:
        writer.writerow(row)

print(f"Zapisano {len(filtered_data)} rekordów do mov_results.csv")
