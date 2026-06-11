import os
import csv
import statistics

os.makedirs("/results", exist_ok=True)

rows = [
    ["CD8T_ISG_1", "ISG15_positive_CD8_T", 8.2, 7.9, 7.5, 7.8, 5.2, 4.9, 3.2, 2.8],
    ["CD8T_ISG_2", "ISG15_positive_CD8_T", 7.8, 7.4, 7.0, 7.6, 5.0, 4.6, 3.0, 2.6],
    ["CD8T_EFF_1", "Effector_CD8_T", 1.2, 1.1, 1.0, 7.5, 7.8, 7.4, 2.0, 1.8],
    ["CD8T_TEX_1", "Exhausted_CD8_T", 3.0, 2.8, 2.7, 6.8, 3.1, 2.8, 7.5, 7.8],
]

out = []
for r in rows:
    cell, cell_type = r[0], r[1]
    isg_score = statistics.mean(r[2:5])
    cytotoxic_score = statistics.mean(r[5:8])
    exhaustion_score = statistics.mean(r[8:10])
    status = "ISG15_high" if isg_score >= 5 else "ISG15_low"
    out.append([cell, cell_type, round(isg_score, 3), round(cytotoxic_score, 3), round(exhaustion_score, 3), status])

with open("/results/toy_ISG15_T_cell_scores.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["cell", "cell_type", "ISG15_signature_score", "cytotoxic_score", "exhaustion_score", "ISG15_status"])
    writer.writerows(out)

with open("/results/toy_ISG15_T_cell_summary.txt", "w") as f:
    f.write("Toy analysis for ISG15-positive T cell signature scoring\n")
    f.write("This run demonstrates ISG15, cytotoxic, and exhaustion signature scoring using a toy dataset.\n")

print("Toy ISG15 T cell analysis completed successfully.")
