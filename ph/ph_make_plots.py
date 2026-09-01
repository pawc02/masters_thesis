from sage.all import *
import csv
import matplotlib.pyplot as plt

# =====================================
# WCZYTANIE DANYCH
# =====================================

data = []

with open("ph_results.csv", newline="") as f:

    reader = csv.DictReader(f)

    for row in reader:

        data.append({
            "smooth_bits": int(row["smooth_bits"]),
            "smooth_time": float(row["smooth_time"]),
            "prime_time": float(row["prime_time"]),
            "speedup": float(row["speedup"])
        })

print("Wczytano", len(data), "rekordów")

# =====================================
# DANE
# =====================================

bits_vals = [d["smooth_bits"] for d in data]

smooth_vals = [d["smooth_time"] for d in data]

prime_vals = [d["prime_time"] for d in data]

speedup_vals = [d["speedup"] for d in data]

# =====================================
# WYKRES CZASU
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    bits_vals,
    smooth_vals,
    s=40,
    alpha=0.7,
    label="Rząd gładki"
)

plt.scatter(
    bits_vals,
    prime_vals,
    s=40,
    alpha=0.7,
    label="Rząd pierwszy"
)

plt.xlabel("Rozmiar podgrupy [bity]", fontsize=14)
plt.ylabel("Czas wykonania [s]", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.legend(fontsize=18)

plt.tight_layout()

plt.savefig(
    "ph_time.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =====================================
# WYKRES PRZYSPIESZENIA
# =====================================

plt.figure(figsize=(8,6))

plt.scatter(
    bits_vals,
    speedup_vals,
    s=40,
    alpha=0.7
)

plt.xlabel("Rozmiar podgrupy [bity]", fontsize=14)
plt.ylabel("Przyspieszenie", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(
    "ph_speedup.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =====================================
# PH_SUMMARY.TXT
# =====================================

with open("ph_summary.txt", "w") as f:

    f.write("=== DATASET ===\n")
    f.write(f"Number of curves = {len(data)}\n\n")

    f.write(f"min(bits) = {min(bits_vals)}\n")
    f.write(f"max(bits) = {max(bits_vals)}\n\n")

    f.write(f"min(smooth time) = {min(smooth_vals)}\n")
    f.write(f"max(smooth time) = {max(smooth_vals)}\n")

    f.write(f"min(prime time) = {min(prime_vals)}\n")
    f.write(f"max(prime time) = {max(prime_vals)}\n\n")

    f.write(f"min speedup = {min(speedup_vals)}\n")
    f.write(f"max speedup = {max(speedup_vals)}\n")
    