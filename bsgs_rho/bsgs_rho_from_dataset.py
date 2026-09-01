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

results = []

bsgs_ok = 0
rho_ok = 0

with open("bsgs_rho_dataset.csv") as f:

    reader = csv.DictReader(f)

    for row in reader:

        p = Integer(row["p"])
        a = Integer(row["a"])
        b = Integer(row["b"])
        N = Integer(row["N"])
        m = Integer(row["m"])

        print("\n================================")
        print("m =", m)
        print("p =", p)

        F = GF(p)
        E = EllipticCurve(F, [a, b])

        # punkt rzędu m
        P = (N // m) * E.random_point()

        while P.is_zero() or P.order() != m:
            P = (N // m) * E.random_point()

        # tajny logarytm
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
        # SPRAWDZENIE POPRAWNOŚCI
        # =====================================

        if x_bsgs == x:
            bsgs_ok += 1

        if x_rho == x:
            rho_ok += 1

        # =====================================
        # ZAPIS
        # =====================================

        results.append({
            "m": int(m),
            "m_bits": int(row["m_bits"]),
            "p": int(p),
            "p_bits": int(row["p_bits"]),
            "a": int(a),
            "b": int(b),
            "N": int(N),
            "bsgs_time": float(bsgs_time),
            "rho_time": float(rho_time),
            "bsgs_table_size": int(ceil(sqrt(m)))
        })


# =====================================
# ZAPIS DO CSV
# =====================================

results.sort(key=lambda d: d["m"])

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

    for row in results:
        writer.writerow(row)

print(f"Współczynnik sukcesu BSGS: {bsgs_ok}/{len(results)}")
print(f"Współczynnik sukcesu Pollard rho: {rho_ok}/{len(results)}")
print(f"Zapisano {len(results)} rekordów do bsgs_rho_results.csv")
