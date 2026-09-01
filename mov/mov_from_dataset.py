from sage.all import *
import time
import csv

# =====================================
# EKSPERYMENT
# =====================================

results = []

before_mov_ok = 0
after_mov_ok = 0

with open("mov_dataset.csv") as f:

    reader = csv.DictReader(f)

    for row in reader:

        p = Integer(row["p"])
        N = Integer(row["N"])
        m = Integer(row["m"])

        print("\n================================")
        print("m =", m)
        print("p =", p)

        F = GF(p)
        E = EllipticCurve(F,[1,0])

        # punkt rzędu m
        P = (N // m) * E.random_point()

        while P.is_zero() or P.order() != m:
            P = (N // m) * E.random_point()

        # tajny logarytm
        x = ZZ.random_element(1,m)

        R = x * P

        # =====================================
        # PRZED REDUKCJĄ MOV
        # =====================================

        start = time.time()
        x_before_mov = discrete_log(R, P, operation='+')
        time_before_mov = time.time() - start

        if x_before_mov == x:
            before_mov_ok += 1

        # =====================================
        # PO REDUKCJI MOV
        # =====================================

        k = Mod(p,m).multiplicative_order()

        Ek = E.base_extend(GF(p**k))

        Pk = Ek(P)
        Rk = Ek(R)

        while True:
            Q = Ek.random_point()
            alpha = Pk.tate_pairing(Q, m, k)
            if alpha.multiplicative_order() == m:
                break

        beta = Rk.tate_pairing(Q, m, k)

        start = time.time()
        x_after_mov = discrete_log(beta, alpha)
        time_after_mov = time.time() - start

        if x_after_mov == x:
            after_mov_ok += 1

        # =====================================
        # ZAPIS
        # =====================================

        results.append({
            "m": int(m),
            "m_bits": int(row["m_bits"]),
            "p": int(p),
            "p_bits": int(row["p_bits"]),
            "N": int(N),
            "embedding_degree": int(k),
            "time_before_mov": float(time_before_mov),
            "time_after_mov": float(time_after_mov)
        })


# =====================================
# ZAPIS DO CSV
# =====================================

results.sort(key=lambda d: d["m"])

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
            "N",
            "embedding_degree",
            "time_before_mov",
            "time_after_mov"
        ]
    )

    writer.writeheader()

    for row in results:
        writer.writerow(row)

print(f"Współczynnik sukcesu przed MOV: {before_mov_ok}/{len(results)}")
print(f"Współczynnik sukcesu po MOV: {after_mov_ok}/{len(results)}")
print(f"Zapisano {len(results)} rekordów do mov_results.csv")
