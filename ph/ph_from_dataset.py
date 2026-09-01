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


# ============================================================
# POHLIG–HELLMAN modulo p^e
# ============================================================

def pohlig_hellman_prime_power(P, R, n, p, e):

    x = 0

    # punkt rzędu p
    P0 = (n // p) * P

    for j in range(e):

        # odejmujemy już znalezione cyfry
        temp = R - x * P

        # projekcja
        Qj = (n // (p**(j + 1))) * temp

        # rozwiązujemy DLP w grupie rzędu p
        zj = bsgs_ecdlp(P0, Qj, p)

        if zj is None:
            return None

        # dodajemy kolejną cyfrę
        x += Integer(zj) * (p**j)

    return x


# ============================================================
# POHLIG–HELLMAN dla ECDLP
# ============================================================

def pohlig_hellman_ecdlp(P, R, n):

    # rozkład rzędu grupy
    factors = factor(n)

    residues = []
    moduli = []

    # ========================================================
    # PODPROBLEMY
    # ========================================================

    for p, e in factors:

        pe = p**e

        # rozwiązanie modulo p^e
        x_pe = pohlig_hellman_prime_power(P, R, n, p, e)

        if x_pe is None:
            return None

        residues.append(x_pe)
        moduli.append(pe)

    # ========================================================
    # CRT
    # ========================================================

    x = CRT_list(residues, moduli)

    return Integer(x % n)


# =====================================
# POJEDYNCZY EKSPERYMENT
# =====================================

def run_ph_experiment(P, n):

    x = ZZ.random_element(1, n)

    R = x * P

    start = time.time()

    x_rec = pohlig_hellman_ecdlp(P, R, n)

    elapsed = time.time() - start

    return elapsed, (x_rec * P == R)

# =====================================
# EKSPERYMENT
# =====================================

results = []

smooth_ok = 0
prime_ok = 0

with open("ph_dataset.csv") as f:

    reader = csv.DictReader(f)

    for row in reader:

        print("\n================================")

        # =====================================
        # SMOOTH CASE
        # =====================================

        p_s = Integer(row["smooth_p"])
        a_s = Integer(row["smooth_a"])
        b_s = Integer(row["smooth_b"])

        N_s = Integer(row["smooth_N"])
        n_s = Integer(row["smooth_order"])

        F_s = GF(p_s)
        E_s = EllipticCurve(F_s, [a_s, b_s])

        P_s = (N_s // n_s) * E_s.random_point()

        while P_s.is_zero() or P_s.order() != n_s:
            P_s = (N_s // n_s) * E_s.random_point()

        times_s = []

        for _ in range(5):

            t, ok = run_ph_experiment(P_s, n_s)

            if ok:
                smooth_ok += 1

            times_s.append(t)

        t_s = sum(times_s) / len(times_s)

        # =====================================
        # PRIME CASE
        # =====================================

        p_p = Integer(row["prime_p"])
        a_p = Integer(row["prime_a"])
        b_p = Integer(row["prime_b"])

        N_p = Integer(row["prime_N"])
        n_p = Integer(row["prime_order"])

        F_p = GF(p_p)
        E_p = EllipticCurve(F_p, [a_p, b_p])

        P_p = (N_p // n_p) * E_p.random_point()

        while P_p.is_zero() or P_p.order() != n_p:
            P_p = (N_p // n_p) * E_p.random_point()

        times_p = []

        for _ in range(5):

            t, ok = run_ph_experiment(P_p, n_p)

            if ok:
                prime_ok += 1

            times_p.append(t)

        t_p = sum(times_p) / len(times_p)

        # =====================================
        # ZAPIS
        # =====================================

        print("smooth bits =", row["smooth_bits"])
        print("smooth time =", t_s)
        print("prime time =", t_p)
        print("speedup =", t_p / t_s)

        results.append({
            "smooth_p": int(p_s),
            "smooth_a": int(a_s),
            "smooth_b": int(b_s),
            "smooth_N": int(N_s),
            "smooth_order": int(n_s),
            "smooth_bits": int(row["smooth_bits"]),
            "smooth_largest_factor": int(row["smooth_largest_factor"]),
            "smooth_largest_factor_bits": int(row["smooth_largest_factor_bits"]),
            "smooth_num_factors": int(row["smooth_num_factors"]),
            "smooth_time": float(t_s),
            "prime_p": int(p_p),
            "prime_a": int(a_p),
            "prime_b": int(b_p),
            "prime_N": int(N_p),
            "prime_order": int(n_p),
            "prime_bits": int(row["prime_bits"]),
            "prime_time": float(t_p),
            "speedup": float(t_p / t_s)
        })


# =====================================
# ZAPIS DO CSV
# =====================================

results.sort(key=lambda d: d["smooth_order"])

with open(
    "ph_results.csv",
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "smooth_p",
            "smooth_a",
            "smooth_b",
            "smooth_N",
            "smooth_order",
            "smooth_bits",
            "smooth_largest_factor",
            "smooth_largest_factor_bits",
            "smooth_num_factors",
            "smooth_time",
            "prime_p",
            "prime_a",
            "prime_b",
            "prime_N",
            "prime_order",
            "prime_bits",
            "prime_time",
            "speedup"
        ]
    )

    writer.writeheader()

    for row in results:
        writer.writerow(row)

print(f"Współczynnik sukcesu Smooth: {smooth_ok}/{5*len(results)}")
print(f"Współczynnik sukcesu Prime: {prime_ok}/{5*len(results)}")
print(f"Zapisano {len(results)} rekordów do ph_results.csv")
