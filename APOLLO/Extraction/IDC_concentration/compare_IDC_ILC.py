import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu, shapiro, normaltest

# Load data
idc = pd.read_csv('Concentration_summary_IDC.csv')
ilc = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Extraction/Concentration_summary_4.12.csv')

def clean_conc(df, col):
    return np.log10([100.0 if str(x).strip() == '<100' else float(x) for x in df[col]])

# Clean concentrations with <100 as 100 for IDC
idc_conc = [100.0 if str(x).strip() == '<100' else float(x) for x in idc['Concentration']]
ilc_conc = [100.0 if str(x).strip() == '<100' else float(x) for x in ilc['Concentration(ng/mL)']]

# Calculate log10 values
idc_log = np.log10(idc_conc)
ilc_log = np.log10(ilc_conc)

# Calculate original scale means and medians
idc_mean_orig = np.mean(idc_conc)
ilc_mean_orig = np.mean(ilc_conc)
idc_median_orig = np.median(idc_conc)
ilc_median_orig = np.median(ilc_conc)

data = [idc_log, ilc_log]
labels = ['IDC', 'ILC']

# Violin plot with dots
plt.figure(figsize=(10, 8))
parts = plt.violinplot(data, showmedians=True, showmeans=True)
plt.xticks([1,2], labels)
plt.ylabel('Log10 Concentration (ng/mL)')
plt.title('Violin plot of IDC vs ILC ctDNA Concentration (log10 scale)')

# Add dots for individual samples
for i, group in enumerate(data):
    x = np.random.normal(i+1, 0.08, size=len(group))  # jitter
    plt.scatter(x, group, color='black', alpha=0.5, s=10)

# t-test
t_stat, t_p = ttest_ind(idc_log, ilc_log, equal_var=False)
# Mann-Whitney U test (Wilcoxon rank-sum)
u_stat, u_p = mannwhitneyu(idc_log, ilc_log, alternative='two-sided')

# Calculate y-axis limits with more space for annotation
y_min = min(min(idc_log), min(ilc_log))
y_max = max(max(idc_log), max(ilc_log))
y_range = y_max - y_min

# Annotate p-value and means below x-axis with more space
plt.text(1.5, y_min - y_range * 0.2, 
         f'Mann-Whitney U test p = {u_p:.3g}', 
         ha='center', color='black', fontsize=12)

# Adjust y-axis limits to accommodate the annotation
plt.ylim(y_min - y_range * 0.25, y_max + y_range * 0.05)

# Adjust layout with more space at the bottom
plt.subplots_adjust(bottom=0.25)

plt.savefig('violin_IDC_vs_ILC_log10.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot histogram of log10 concentrations
plt.figure(figsize=(8,6))
plt.hist(idc_log, bins=30, alpha=0.6, label='IDC', color='skyblue', edgecolor='black')
plt.hist(ilc_log, bins=30, alpha=0.6, label='ILC', color='orange', edgecolor='black')
plt.xlabel('Log10 Concentration (ng/mL)')
plt.ylabel('Frequency')
plt.title('Histogram of Log10 Concentration: IDC vs ILC')
plt.legend()
plt.tight_layout()
plt.savefig('histogram_IDC_vs_ILC_log10.png', dpi=300)
plt.close()

# Normality tests
shapiro_idc = shapiro(idc_log)
shapiro_ilc = shapiro(ilc_log)
normtest_idc = normaltest(idc_log)
normtest_ilc = normaltest(ilc_log)

# Summary statistics
def summary_stats(arr):
    return np.median(arr), np.mean(arr), np.std(arr)

idc_median, idc_mean, idc_sd = summary_stats(idc_log)
ilc_median, ilc_mean, ilc_sd = summary_stats(ilc_log)

# Cohen's d
def cohens_d(x, y):
    nx, ny = len(x), len(y)
    pooled_sd = np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / (nx+ny-2))
    return (np.mean(x) - np.mean(y)) / pooled_sd

d = cohens_d(idc_log, ilc_log)

# Print results
print(f'T-test: t = {t_stat:.3f}, p = {t_p:.3g}')
print(f'Mann-Whitney U test: U = {u_stat:.3f}, p = {u_p:.3g}')
print(f'IDC: median = {idc_median:.3f}, mean = {idc_mean:.3f}, sd = {idc_sd:.3f}')
print(f'ILC: median = {ilc_median:.3f}, mean = {ilc_mean:.3f}, sd = {ilc_sd:.3f}')
print(f"Cohen's d effect size: {d:.3f}")
print(f'Shapiro-Wilk normality test (IDC): W = {shapiro_idc.statistic:.3f}, p = {shapiro_idc.pvalue:.3g}')
print(f'Shapiro-Wilk normality test (ILC): W = {shapiro_ilc.statistic:.3f}, p = {shapiro_ilc.pvalue:.3g}')
print(f"D'Agostino K^2 normality test (IDC): stat = {normtest_idc.statistic:.3f}, p = {normtest_idc.pvalue:.3g}")
print(f"D'Agostino K^2 normality test (ILC): stat = {normtest_ilc.statistic:.3f}, p = {normtest_ilc.pvalue:.3g}")

# Mann-Whitney U test (Wilcoxon rank-sum) on original concentrations
u_stat_orig, u_p_orig = mannwhitneyu(idc_conc, ilc_conc, alternative='two-sided')
print(f'Mann-Whitney U test (original concentration): U = {u_stat_orig:.3f}, p = {u_p_orig:.3g}')

# Add annotations below x-axis with mean and median values for both cohorts
plt.figure(figsize=(8,6))
plt.bar(0, idc_mean, color='skyblue', label='IDC')
plt.bar(1, ilc_mean, color='orange', label='ILC')
plt.xticks([0, 1], ['IDC', 'ILC'])
plt.ylabel('Concentration (ng/mL)')
plt.title('Mean Concentration: IDC vs ILC')
plt.legend()
plt.tight_layout()
plt.savefig('bar_IDC_vs_ILC_mean.png', dpi=300)
plt.close()

plt.figure(figsize=(8,6))
plt.bar(0, idc_median, color='skyblue', label='IDC')
plt.bar(1, ilc_median, color='orange', label='ILC')
plt.xticks([0, 1], ['IDC', 'ILC'])
plt.ylabel('Concentration (ng/mL)')
plt.title('Median Concentration: IDC vs ILC')
plt.legend()
plt.tight_layout()
plt.savefig('bar_IDC_vs_ILC_median.png', dpi=300)
plt.close() 