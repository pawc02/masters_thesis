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


# ============================================================
# LOSOWA KRZYWA
# ============================================================

def random_curve(bits):

    p = random_prime(2**bits - 1,
                     lbound=2**(bits-1))

    F = GF(p)

    while True:

        a = F.random_element()
        b = F.random_element()

        if 4*a**3 + 27*b**2 != 0:
            break

    E = EllipticCurve(F, [a, b])

    return E


# ============================================================
# FIND SMOOTH SUBGROUP
# ============================================================

def find_smooth_subgroup(bits,
                         smooth_bound=10000,
                         target_bits=30):

    while True:

        E = random_curve(bits)

        N = E.cardinality()

        fac = factor(N)

        smooth_part = Integer(1)

        for p, e in fac:

            if p < smooth_bound:

                smooth_part *= p**e

        if smooth_part.nbits() != target_bits:
            continue

        print("\nFOUND SMOOTH GROUP")
        print("N =", N)
        print("smooth part =", smooth_part)
        print("bitlength =", smooth_part.nbits())

        largest_factor = max([p for p,e in factor(smooth_part)])

        print("largest prime factor =", largest_factor)
        print("largest factor bits =", largest_factor.nbits())

        for _ in range(1000):

            Q = E.random_point()

            P = (N // smooth_part) * Q

            if P.order() == smooth_part:
                return E, P, smooth_part

        # nie znaleziono punktu rzędu smooth_part
        # wracamy do początku while True


# ============================================================
# FIND PRIME SUBGROUP
# ============================================================

def find_prime_subgroup(bits,
                        target_bits):

    while True:

        E = random_curve(bits)

        N = E.cardinality()

        fac = factor(N)

        primes = [p for p,e in fac if p.is_prime()]

        if len(primes) == 0:
            continue

        q = max(primes)

        if q.nbits() != target_bits:
            continue

        print("\nFOUND PRIME GROUP")
        print("N =", N)
        print("prime =", q)
        print("bitlength =", q.nbits())
        print("largest prime factor =", q)
        
        while True:

            Q = E.random_point()

            P = (N // q) * Q

            if P.order() == q:
                return E, P, q


# ============================================================
# POJEDYNCZY EKSPERYMENT
# ============================================================

def run_ph_experiment(P, n):

    x = ZZ.random_element(1, n)

    R = x * P

    start = time.time()

    x_rec = pohlig_hellman_ecdlp(P, R, n)

    elapsed = time.time() - start

    ok = (x_rec * P == R)

    return elapsed, ok


# ============================================================
# EKSPERYMENT
# ============================================================

def experiment(bits):

    print("\n================================")
    print("TARGET BITS =", bits)
    print("================================")

    # =====================================
    # SMOOTH CASE
    # =====================================

    E_s, P_s, n_s = find_smooth_subgroup(
        bits=bits,
        smooth_bound=10000,
        target_bits=bits
    )

    times_s = []

    for _ in range(5):

        t, ok = run_ph_experiment(P_s, n_s)

        if not ok:
            raise ValueError("PH failed (smooth)")

        times_s.append(t)

    t_s = sum(times_s) / len(times_s)

    p_s = E_s.base_field().order()
    a_s = E_s.a4()
    b_s = E_s.a6()
    N_s = E_s.cardinality()
    largest_factor = max([p for p,e in factor(n_s)])

    # =====================================
    # PRIME CASE
    # =====================================

    E_p, P_p, n_p = find_prime_subgroup(
        bits=bits,
        target_bits=n_s.nbits()
    )

    times_p = []

    for _ in range(5):

        t, ok = run_ph_experiment(P_p, n_p)

        if not ok:
            raise ValueError("PH failed (prime)")

        times_p.append(t)

    t_p = sum(times_p) / len(times_p)

    print("smooth time =", t_s)
    print("prime time =", t_p)

    p_p = E_p.base_field().order()
    a_p = E_p.a4()
    b_p = E_p.a6()
    N_p = E_p.cardinality()

    return {
        "smooth_p": int(p_s),
        "smooth_a": int(a_s),
        "smooth_b": int(b_s),
        "smooth_N": int(N_s),
        "smooth_order": int(n_s),
        "smooth_bits": int(n_s.nbits()),
        "smooth_largest_factor": int(largest_factor),
        "smooth_largest_factor_bits": int(largest_factor.nbits()),
        "smooth_num_factors": len(factor(n_s)),
        "smooth_time": float(t_s),
        "prime_p": int(p_p),
        "prime_a": int(a_p),
        "prime_b": int(b_p),
        "prime_N": int(N_p),
        "prime_order": int(n_p),
        "prime_bits": int(n_p.nbits()),
        "prime_time": float(t_p),
        "speedup": float(t_p / t_s)
    }


# =====================================
# ZBIERANIE DANYCH
# =====================================

all_data = []

for bits in range(26,36):

    for _ in range(5):

        try:

            res = experiment(bits)

            all_data.append(res)

        except Exception as e:

            print("ERROR:", e)

all_data.sort(key=lambda d: d["smooth_order"])

# =====================================
# DATASET
# =====================================

with open(
    "ph_dataset.csv",
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
            "prime_p",
            "prime_a",
            "prime_b",
            "prime_N",
            "prime_order",
            "prime_bits"
        ]
    )

    writer.writeheader()

    for row in all_data:

        writer.writerow({
            "smooth_p": row["smooth_p"],
            "smooth_a": row["smooth_a"],
            "smooth_b": row["smooth_b"],
            "smooth_N": row["smooth_N"],
            "smooth_order": row["smooth_order"],
            "smooth_bits": row["smooth_bits"],
            "smooth_largest_factor": row["smooth_largest_factor"],
            "smooth_largest_factor_bits": row["smooth_largest_factor_bits"],
            "smooth_num_factors": row["smooth_num_factors"],
            "prime_p": row["prime_p"],
            "prime_a": row["prime_a"],
            "prime_b": row["prime_b"],
            "prime_N": row["prime_N"],
            "prime_order": row["prime_order"],
            "prime_bits": row["prime_bits"]
        })

print(f"Zapisano {len(all_data)} rekordów do ph_dataset.csv")

# =====================================
# WYNIKI EKSPERYMENTÓW
# =====================================

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

    for row in all_data:
        writer.writerow(row)

print(f"Zapisano {len(all_data)} rekordów do ph_results.csv")
