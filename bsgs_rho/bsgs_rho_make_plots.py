from sage.all import *
import csv
import matplotlib.pyplot as plt
import math

# =====================================
# WCZYTANIE DANYCH
# =====================================

data = []

with open("bsgs_rho_results.csv", newline="") as f:

    reader = csv.DictReader(f)

    for row in reader:

        data.append({

            "m": Integer(row["m"]),
            "m_bits": Integer(row["m_bits"]),

            "bsgs_time": float(row["bsgs_time"]),
            "rho_time": float(row["rho_time"]),

            "bsgs_table_size":
                Integer(row["bsgs_table_size"])
        })

print("Wczytano", len(data), "rekordów")

data.sort(key=lambda d: d["m"])

# =====================================
# DANE
# =====================================

m_vals = [d["m"] for d in data]

bsgs_vals = [d["bsgs_time"] for d in data]
rho_vals = [d["rho_time"] for d in data]

table_vals = [d["bsgs_table_size"] for d in data]

log_m_vals = [math.log2(float(m)) for m in m_vals]

# =====================================
# REGRESJA LINIOWA
# =====================================

def linear_regression(points):

    X = matrix([[x, 1] for x, y in points])
    Y = vector([y for x, y in points])

    sol = ((X.transpose() * X).inverse() * X.transpose() * Y)

    return sol[0], sol[1]

# =====================================
# WYKRES CZASU
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    log_m_vals,
    bsgs_vals,
    s=40,
    alpha=0.7,
    label="BSGS"
)

plt.scatter(
    log_m_vals,
    rho_vals,
    s=40,
    alpha=0.7,
    label="rho Pollarda"
)

plt.xlabel("Rozmiar podgrupy [bity]", fontsize=14)
plt.ylabel("Czas wykonania [s]", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.legend(fontsize=18)

plt.tight_layout()

plt.savefig("bsgs_rho_time.png", dpi=300, bbox_inches="tight")

plt.close()

# =====================================
# LOG-LOG
# =====================================

points_bsgs = [
    (math.log2(float(m)), math.log2(float(t)))
    for m, t in zip(m_vals, bsgs_vals)
]

points_rho = [
    (math.log2(float(m)), math.log2(float(t)))
    for m, t in zip(m_vals, rho_vals)
]

# =====================================
# REGRESJA
# =====================================

a_bsgs, b_bsgs = linear_regression(points_bsgs)
a_rho, b_rho = linear_regression(points_rho)

print("\n=== LOG-LOG FIT ===")

print("BSGS slope =", float(a_bsgs))
print("Pollard rho slope =", float(a_rho))

# =====================================
# LINIE DOPASOWANIA
# =====================================

fit_bsgs = [
    (x, a_bsgs*x + b_bsgs)
    for x,_ in points_bsgs
]

fit_rho = [
    (x, a_rho*x + b_rho)
    for x,_ in points_rho
]

# =====================================
# WYKRES LOG-LOG
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    [x for x,y in points_bsgs],
    [y for x,y in points_bsgs],
    s=40,
    alpha=0.7,
    label="BSGS"
)

plt.scatter(
    [x for x,y in points_rho],
    [y for x,y in points_rho],
    s=40,
    alpha=0.7,
    label="rho Pollarda"
)

plt.plot(
    [x for x,y in fit_bsgs],
    [y for x,y in fit_bsgs],
    linewidth=2
)

plt.plot(
    [x for x,y in fit_rho],
    [y for x,y in fit_rho],
    linewidth=2
)

plt.xlabel("log₂(m)", fontsize=14)
plt.ylabel("log₂(czas wykonania)", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.legend(fontsize=18)

plt.tight_layout()

plt.savefig(
    "bsgs_rho_loglog.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =====================================
# T(m)/sqrt(m)
# =====================================

ratios_bsgs = [
    t / math.sqrt(float(m))
    for m, t in zip(m_vals, bsgs_vals)
]

ratios_rho = [
    t / math.sqrt(float(m))
    for m, t in zip(m_vals, rho_vals)
]

avg_bsgs = sum(ratios_bsgs) / len(ratios_bsgs)
avg_rho  = sum(ratios_rho) / len(ratios_rho)

print("\n=== T(m)/sqrt(m) ===")

print("avg BSGS =", float(avg_bsgs))
print("avg rho  =", float(avg_rho))

print("ratio BSGS/rho =", float(avg_bsgs / avg_rho))

# =====================================
# BSGS_RHO_SUMMARY.TXT
# =====================================

with open("bsgs_rho_summary.txt", "w") as f:

    f.write("=== DATASET ===\n")
    f.write(f"Number of curves = {len(data)}\n\n")
    f.write(f"min(m) = {min(m_vals)}\n")
    f.write(f"max(m) = {max(m_vals)}\n")
    f.write(f"min(log2(m)) = {float(min(log_m_vals))}\n")
    f.write(f"max(log2(m)) = {float(max(log_m_vals))}\n\n")
    f.write(f"min(BSGS time) = {min(bsgs_vals)}\n")
    f.write(f"max(BSGS time) = {max(bsgs_vals)}\n")
    f.write(f"min(rho time) = {min(rho_vals)}\n")
    f.write(f"max(rho time) = {max(rho_vals)}\n\n")
    f.write("=== LOG-LOG FIT ===\n")
    f.write(f"BSGS slope = {float(a_bsgs)}\n")
    f.write(f"Pollard rho slope = {float(a_rho)}\n\n")
    f.write("=== T(m)/sqrt(m) ===\n")
    f.write(f"avg BSGS = {float(avg_bsgs)}\n")
    f.write(f"avg rho = {float(avg_rho)}\n")
    f.write(f"ratio BSGS/rho = {float(avg_bsgs / avg_rho)}\n\n")

    # =====================================
    # EKSTRAPOLACJA RHO POLLARDA
    # =====================================

    f.write("\n=== POLLARD RHO EXTRAPOLATION ===\n")
    SECONDS_PER_YEAR = 365 * 24 * 3600
    for bits in [128, 192, 256, 512, 795]:
        T_seconds = avg_rho * 2**(bits/2)
        years = T_seconds / SECONDS_PER_YEAR
        f.write(f"{bits} bits: {years:.3e} years\n")

    f.write("\nNFS-DL (literatura):\n")
    f.write("795-bit DLP ≈ 3200 core-years\n")

# =====================================
# WYKRES STAŁEJ
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    log_m_vals,
    ratios_bsgs,
    s=40,
    alpha=0.7,
    label='BSGS'
)

plt.scatter(
    log_m_vals,
    ratios_rho,
    s=40,
    alpha=0.7,
    label='rho Pollarda'
)

plt.xlabel('Rozmiar podgrupy [bity]', fontsize=14)
plt.ylabel(r'$T(m)/\sqrt{m}$', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.legend(fontsize=18)

plt.tight_layout()

plt.savefig(
    'bsgs_rho_ratio.png',
    dpi=300,
    bbox_inches='tight'
)

plt.close()

# =====================================
# ZAKOŃCZENIE
# =====================================

print("\nSaved plots:")

print(" - bsgs_rho_time.png")
print(" - bsgs_rho_loglog.png")
print(" - bsgs_rho_ratio.png")

print("\nSaved summary:")
print(" - bsgs_rho_summary.txt")
