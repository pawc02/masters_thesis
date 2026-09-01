from sage.all import *
import time
import csv

# =====================================
# BABY STEP GIANT STEP
# =====================================

def bsgs_ecdlp(P, R, n):

    m = ceil(sqrt(n))

    # =====================================
    # BABY STEPS
    # =====================================

    baby_table = {}

    current = P.curve()(0)

    for j in range(m):

        baby_table[current] = j
        current += P

    # =====================================
    # GIANT STEPS
    # =====================================

    giant_step = m * P

    current = R

    for i in range(m):

        if current in baby_table:

            j = baby_table[current]
            return i * m + j

        current -= giant_step

    return None


# =====================================
# RHO POLLARDA
# =====================================

def partition(X):

    return Integer(X[0]) % 3


def rho_step(X, a, b, P, R, m):

    subset = partition(X)

    # S1
    if subset == 0:

        X = X + P
        a = (a + 1) % m

    # S2
    elif subset == 1:

        X = 2 * X
        a = (2 * a) % m
        b = (2 * b) % m

    # S3
    else:

        X = X + R
        b = (b + 1) % m

    return X, a, b


def pollard_rho_ecdlp(P, R, m):

    while True:

        # =====================================
        # LOSOWY START
        # =====================================

        a = ZZ.random_element(m)
        b = ZZ.random_element(m)

        X = a*P + b*R

        # tortoise
        X1, a1, b1 = X, a, b

        # hare
        X2, a2, b2 = X, a, b

        # =====================================
        # FLOYD CYCLE FINDING
        # =====================================

        while True:

            # tortoise
            X1, a1, b1 = rho_step(X1, a1, b1, P, R, m)

            # hare
            X2, a2, b2 = rho_step(X2, a2, b2, P, R, m)
            X2, a2, b2 = rho_step(X2, a2, b2, P, R, m)

            # collision
            if X1 == X2:

                s = (a1 - a2) % m
                t = (b2 - b1) % m

                if t == 0:
                    break

                if gcd(t, m) != 1:
                    break

                try:

                    x = (s * inverse_mod(t, m)) % m
                    return x

                except:
                    break


# =====================================
# EKSPERYMENT
# =====================================

def experiment(bits=32):

    print("\n==== bits =", bits, "====")

    # =====================================
    # LOSOWE p
    # =====================================

    p = random_prime(2**bits - 1, lbound=2**(bits-1))

    F = GF(p)

    # =====================================
    # LOSOWA KRZYWA
    # =====================================

    while True:

        a = F.random_element()
        b = F.random_element()

        if 4*a**3 + 27*b**2 != 0:
            break

    E = EllipticCurve(F, [a, b])

    print("Curve:")
    print(f"y^2 = x^3 + {a}x + {b}")

    # =====================================
    # RZĄD KRZYWEJ
    # =====================================

    N = E.cardinality()

    print("Group order N =", N)

    # =====================================
    # DUŻY CZYNNIK PIERWSZY
    # =====================================

    factors = factor(N)

    m = max([f[0] for f in factors])

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
    # BSGS
    # =====================================

    start = time.time()

    x_bsgs = bsgs_ecdlp(P, R, m)

    bsgs_time = time.time() - start

    # =====================================
    # RHO POLLARDA
    # =====================================

    start = time.time()

    x_rho = pollard_rho_ecdlp(P, R, m)

    rho_time = time.time() - start

    # =====================================
    # WERYFIKACJA
    # =====================================

    print("BSGS time =", bsgs_time)
    print("BSGS correct? =", x_bsgs == x)

    print("Pollard rho time =", rho_time)
    print("Pollard rho correct? =", x_rho == x)

    # =====================================
    # TABLICA BSGS
    # =====================================

    table_size = ceil(sqrt(m))

    print("baby table size =", table_size)

    return {
        "m": int(m),
        "m_bits": int(m.nbits()),
        "p": int(p),
        "p_bits": int(p.nbits()),
        "a": int(a),
        "b": int(b),
        "N": int(N),
        "bsgs_time": float(bsgs_time),
        "rho_time": float(rho_time),
        "bsgs_table_size": int(table_size)
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
    "bsgs_rho_dataset.csv",
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

print(f"Zapisano {len(filtered_data)} rekordów do bsgs_rho_dataset.csv")

# =====================================
# WYNIKI EKSPERYMENTÓW
# =====================================

with open(
    "bsgs_rho_results.csv",
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
            "bsgs_time",
            "rho_time",
            "bsgs_table_size"
        ]
    )

    writer.writeheader()

    for row in filtered_data:
        writer.writerow(row)

print(f"Zapisano {len(filtered_data)} rekordów do bsgs_rho_results.csv")
