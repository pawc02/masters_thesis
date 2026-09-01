from sage.all import *
import csv
import matplotlib.pyplot as plt
import math

# =====================================
# WCZYTANIE DANYCH
# =====================================

data = []

with open("mov_results.csv", newline="") as f:

    reader = csv.DictReader(f)

    for row in reader:

        data.append({

            "m": Integer(row["m"]),
            "m_bits": Integer(row["m_bits"]),

            "time_before_mov": float(row["time_before_mov"]),
            "time_after_mov": float(row["time_after_mov"])
        })

print("Wczytano", len(data), "rekordów")

data.sort(key=lambda d: d["m"])

# =====================================
# DANE
# =====================================

m_vals = [d["m"] for d in data]

before_vals = [d["time_before_mov"] for d in data]
after_vals = [d["time_after_mov"] for d in data]

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
    before_vals,
    s=40,
    alpha=0.7,
    label="Przed MOV"
)

plt.scatter(
    log_m_vals,
    after_vals,
    s=40,
    alpha=0.7,
    label="Po MOV"
)

plt.xlabel("Rozmiar podgrupy [bity]", fontsize=14)
plt.ylabel("Czas wykonania [s]", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.legend(fontsize=18)

plt.tight_layout()

plt.savefig("mov_time.png", dpi=300, bbox_inches="tight")

plt.close()

# =====================================
# LOG-LOG
# =====================================

points_before = [
    (math.log2(float(m)), math.log2(float(t)))
    for m, t in zip(m_vals, before_vals)
]

points_after = [
    (math.log2(float(m)), math.log2(float(t)))
    for m, t in zip(m_vals, after_vals)
]

# =====================================
# REGRESJA
# =====================================

a_before, b_before = linear_regression(points_before)
a_after, b_after = linear_regression(points_after)

print("\n=== LOG-LOG FIT ===")
print("Przed MOV slope =", float(a_before))
print("Po MOV slope =", float(a_after))

# =====================================
# LINIE DOPASOWANIA
# =====================================

fit_before = [
    (x, a_before*x + b_before)
    for x,_ in points_before
]

fit_after = [
    (x, a_after*x + b_after)
    for x,_ in points_after
]

# =====================================
# WYKRES LOG-LOG
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    [x for x,y in points_before],
    [y for x,y in points_before],
    s=40,
    alpha=0.7,
    label="Przed MOV"
)

plt.scatter(
    [x for x,y in points_after],
    [y for x,y in points_after],
    s=40,
    alpha=0.7,
    label="Po MOV"
)

plt.plot(
    [x for x,y in fit_before],
    [y for x,y in fit_before],
    linewidth=2
)

plt.plot(
    [x for x,y in fit_after],
    [y for x,y in fit_after],
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
    "mov_loglog.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =====================================
# T(m)/sqrt(m)
# =====================================

ratios_before = [
    t / math.sqrt(float(m))
    for m, t in zip(m_vals, before_vals)
]

ratios_after = [
    t / math.sqrt(float(m))
    for m, t in zip(m_vals, after_vals)
]

avg_before = sum(ratios_before) / len(ratios_before)
avg_after = sum(ratios_after) / len(ratios_after)


print("\n=== T(m)/sqrt(m) ===")

print("avg before MOV =", float(avg_before))
print("avg after MOV =", float(avg_after))

print("ratio before/after =", float(avg_before / avg_after))

# =====================================
# MOV_SUMMARY.TXT
# =====================================

with open("mov_summary.txt", "w") as f:

    f.write("=== DATASET ===\n")
    f.write(f"Number of curves = {len(data)}\n\n")
    f.write(f"min(m) = {min(m_vals)}\n")
    f.write(f"max(m) = {max(m_vals)}\n")
    f.write(f"min(before MOV time) = {min(before_vals)}\n")
    f.write(f"max(before MOV time) = {max(before_vals)}\n")
    f.write(f"min(after MOV time) = {min(after_vals)}\n")
    f.write(f"max(after MOV time) = {max(after_vals)}\n\n")
    f.write("=== LOG-LOG FIT ===\n")
    f.write(f"before MOV slope = {float(a_before)}\n")
    f.write(f"after MOV slope = {float(a_after)}\n\n")
    f.write("=== T(m)/sqrt(m) ===\n")
    f.write(f"avg before MOV = {float(avg_before)}\n")
    f.write(f"avg after MOV = {float(avg_after)}\n")
    f.write(f"ratio before/after = {float(avg_before/avg_after)}\n")

# =====================================
# WYKRES STAŁEJ
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    log_m_vals,
    ratios_before,
    s=40,
    alpha=0.7,
    label="Przed MOV"
)

plt.scatter(
    log_m_vals,
    ratios_after,
    s=40,
    alpha=0.7,
    label="Po MOV"
)

plt.xlabel("Rozmiar podgrupy [bity]", fontsize=14)
plt.ylabel(r"$T(m)/\sqrt{m}$", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.legend(fontsize=18)

plt.tight_layout()

plt.savefig(
    "mov_ratio.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =====================================
# ZAKOŃCZENIE
# =====================================

print("\nSaved plots:")

print(" - mov_time.png")
print(" - mov_loglog.png")
print(" - mov_ratio.png")

print("\nSaved summary:")
print(" - mov_summary.txt")
