# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import numpy as np

# Read the cohort files
idc_df = pd.read_csv('idc_cohort.csv')
ilc_df = pd.read_csv('ilc_cohort.csv')
genomic_file = 'Guardant Project - Genomic Data(UPMC data - Initial Test Only).csv'
genomic_df = pd.read_csv(genomic_file, encoding='latin-1')

# Get column names and IDs
gh_id_col = 'Can delete duplicates of GH_ID'
genomic_gh_id_col = 'ï»¿GH_ID'

# Get sets of patient IDs
idc_ids = set(idc_df[gh_id_col])
ilc_ids = set(ilc_df[gh_id_col])
histology_ids = idc_ids.union(ilc_ids)
genomic_ids = set(genomic_df[genomic_gh_id_col])

# Calculate overlaps
idc_genomic_overlap = idc_ids.intersection(genomic_ids)
ilc_genomic_overlap = ilc_ids.intersection(genomic_ids)

# Create figure with adjusted size and spacing
plt.figure(figsize=(16, 6))
plt.subplots_adjust(wspace=0.3, bottom=0.2)

# 1. Create main Venn diagram
plt.subplot(1, 3, 1)
v = venn2([histology_ids, genomic_ids], 
        set_labels=('Histology Cohort\n(n={})'.format(len(histology_ids)), 
                   'Genomic Cohort\n(n={})'.format(len(genomic_ids))),
        set_colors=('#3498db', '#e67e22'),
        alpha=0.7)

# Add proportion annotations
# For histology cohort
hist_with_genomic = len(histology_ids.intersection(genomic_ids)) / len(histology_ids) * 100
plt.annotate('{:.1f}% of Histology\nhas Genomic Data'.format(hist_with_genomic),
             xy=(-0.5, 0.3), xytext=(-1.0, 0.5),
             arrowprops=dict(facecolor='black', shrink=0.05),
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

# For genomic cohort
genomic_with_hist = len(histology_ids.intersection(genomic_ids)) / len(genomic_ids) * 100
plt.annotate('{:.1f}% of Genomic\nhas Histology Data'.format(genomic_with_hist),
             xy=(0.5, 0.3), xytext=(1.0, 0.5),
             arrowprops=dict(facecolor='black', shrink=0.05),
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

# Add IDC/ILC text in overlap region
overlap_x = 0  # Center of overlap region
overlap_y = -0.4  # Below center of overlap region
plt.text(overlap_x, overlap_y, 
         'IDC: {} ({:.1f}%)\nILC: {} ({:.1f}%)'.format(
             len(idc_genomic_overlap),
             len(idc_genomic_overlap)/len(histology_ids.intersection(genomic_ids))*100,
             len(ilc_genomic_overlap),
             len(ilc_genomic_overlap)/len(histology_ids.intersection(genomic_ids))*100),
         ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8))

plt.title('Patient Overlap Between\nHistology and Genomic Cohorts', pad=20)

# 2. Create IDC/ILC proportion circle for all histology
plt.subplot(1, 3, 2)
plt.pie([len(idc_ids), len(ilc_ids)], 
        labels=['IDC\n(n={})\n{:.1f}%'.format(len(idc_ids), len(idc_ids)/len(histology_ids)*100), 
                'ILC\n(n={})\n{:.1f}%'.format(len(ilc_ids), len(ilc_ids)/len(histology_ids)*100)],
        colors=['#3498db', '#e74c3c'],
        autopct='',
        startangle=90)
plt.title('IDC vs ILC Distribution\nin All Histology Cohort')

# 3. Create IDC/ILC proportion circle for overlapping cohort
plt.subplot(1, 3, 3)
plt.pie([len(idc_genomic_overlap), len(ilc_genomic_overlap)], 
        labels=['IDC\n(n={})\n{:.1f}%'.format(
                    len(idc_genomic_overlap),
                    len(idc_genomic_overlap)/len(histology_ids.intersection(genomic_ids))*100), 
                'ILC\n(n={})\n{:.1f}%'.format(
                    len(ilc_genomic_overlap),
                    len(ilc_genomic_overlap)/len(histology_ids.intersection(genomic_ids))*100)],
        colors=['#3498db', '#e74c3c'],
        autopct='',
        startangle=90)
plt.title('IDC vs ILC Distribution\nin Overlapping Cohort')

# Add numerical summary below the plots in a more compact format
summary_text = (
    "Total Unique Patients: {total} | Histology Only: {hist_only} ({hist_only_pct:.1f}%) | "
    "Genomic Only: {gen_only} ({gen_only_pct:.1f}%) | Both Datasets: {both} ({both_pct:.1f}%)\n"
    "IDC with Genomic: {idc_gen} ({idc_gen_pct:.1f}% of IDC) | "
    "ILC with Genomic: {ilc_gen} ({ilc_gen_pct:.1f}% of ILC)").format(
        total=len(histology_ids.union(genomic_ids)),
        hist_only=len(histology_ids - genomic_ids),
        hist_only_pct=len(histology_ids - genomic_ids)/len(histology_ids.union(genomic_ids))*100,
        gen_only=len(genomic_ids - histology_ids),
        gen_only_pct=len(genomic_ids - histology_ids)/len(histology_ids.union(genomic_ids))*100,
        both=len(histology_ids.intersection(genomic_ids)),
        both_pct=len(histology_ids.intersection(genomic_ids))/len(histology_ids.union(genomic_ids))*100,
        idc_gen=len(idc_genomic_overlap),
        idc_gen_pct=len(idc_genomic_overlap)/len(idc_ids)*100,
        ilc_gen=len(ilc_genomic_overlap),
        ilc_gen_pct=len(ilc_genomic_overlap)/len(ilc_ids)*100)

plt.figtext(0.5, 0.02, summary_text, fontsize=10, ha='center', va='bottom')

# Save with minimal margins
plt.savefig('venn_diagram_compact.png', bbox_inches='tight', dpi=300)
print("\nPlot saved as 'venn_diagram_compact.png'")

# Print numerical summary
print("\nNumerical Summary:")
print("------------------")
print("Total Unique Patients: {}".format(len(histology_ids.union(genomic_ids))))
print("Histology Only: {} ({:.1f}%)".format(
    len(histology_ids - genomic_ids),
    len(histology_ids - genomic_ids)/len(histology_ids.union(genomic_ids))*100))
print("Genomic Only: {} ({:.1f}%)".format(
    len(genomic_ids - histology_ids),
    len(genomic_ids - histology_ids)/len(histology_ids.union(genomic_ids))*100))
print("Both Datasets: {} ({:.1f}%)".format(
    len(histology_ids.intersection(genomic_ids)),
    len(histology_ids.intersection(genomic_ids))/len(histology_ids.union(genomic_ids))*100))
print("  - IDC with Genomic: {} ({:.1f}% of IDC)".format(
    len(idc_genomic_overlap),
    len(idc_genomic_overlap)/len(idc_ids)*100))
print("  - ILC with Genomic: {} ({:.1f}% of ILC)".format(
    len(ilc_genomic_overlap),
    len(ilc_genomic_overlap)/len(ilc_ids)*100))
print("\nCross-Dataset Coverage:")
print("  - {:.1f}% of Histology cohort has Genomic data".format(hist_with_genomic))
print("  - {:.1f}% of Genomic cohort has Histology data".format(genomic_with_hist)) 