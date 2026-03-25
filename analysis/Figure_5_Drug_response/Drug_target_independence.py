import pandas as pd
from datetime import date
import os
import numpy as np
import math
from scipy.stats import norm, fisher_exact
from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns


def odds_ratio_ci(table, odds_ratio: float, alpha : float = 0.05):
    """
    Calculate Odds Ratio and Confidence Interval using the log method.

    table: containing counts about responders and non-responders to a drug target
    alpha: significance level (default 0.05 for 95% CI)
    """

    # If OR is zero or infinity, CI can't be computed normally
    if odds_ratio == 0 or math.isinf(odds_ratio):
        return np.nan, np.nan

    # Convert to float array
    table = np.array(table, dtype=float)

    # Unpack
    a, b = table[0]  # Exposed with/without outcome
    c, d = table[1]  # Not exposed with/without outcome

    # Log OR and Standard Error
    log_or = math.log(odds_ratio)
    SE = math.sqrt(1/a + 1/b + 1/c + 1/d)

    # Z-score for desired CI
    z = norm.ppf(1 - alpha/2)

    # CI for log(OR)
    lower_log = log_or - z * SE
    upper_log = log_or + z * SE

    # Transform back to OR scale
    lower = math.exp(lower_log)
    upper = math.exp(upper_log)

    return lower, upper


def get_stats(df):
    # Compute OR + CI for each endotype vs all others
    # Prepare results DataFrame
    results = []

    # Loop over all pairwise comparisons
    for e1, e2 in combinations(df['Endotype'], 2):
        r1 = df.loc[df['Endotype'] == e1, 'Total_resp'].values[0]
        nr1 = df.loc[df['Endotype'] == e1, 'Total_non'].values[0]
        r2 = df.loc[df['Endotype'] == e2, 'Total_resp'].values[0]
        nr2 = df.loc[df['Endotype'] == e2, 'Total_non'].values[0]

        # Contingency table
        table = [[r1, nr1], [r2, nr2]]

        # Fisher's Exact
        odds_ratio, p_value = fisher_exact(table, alternative='greater')

        # CI
        ci_low, ci_high = odds_ratio_ci(table, odds_ratio)

        # Store results (both directions if needed)
        results.append({
            'Comparison': f'{e1} vs {e2}',
            'Endotype_1': e1,
            'Endotype_2': e2,
            'OR': odds_ratio,
            'CI_low': ci_low,
            'CI_high': ci_high,
            'p_value': p_value
        })

    results_df = pd.DataFrame(results)
    return results_df


def main(save_folder):
    filename = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/raw_data/Table_S9.xlsx'
    # Sample dataframe
    df = pd.DataFrame([
        ['E5',  ['E5', 'E11/E12/E13'],  3, 3, 6, 0, 1, 4], #  'IL17_resp', 'IL23_resp', 'TNFa_resp', 'IL17_non', 'IL23_non', 'TNFa_non'
        ['E11', ['E11', 'E5/E12/E13'],  5, 6, 9, 2, 7, 13],
        ['E12', ['E12', 'E5/E11/E13'],  5, 24, 15, 3, 7, 15],
        ['E13', ['E13', 'E5/E11/E12'],  4, 4, 7, 1, 2, 8],
    ], columns=[
        'Endotype', 'Others',
        'IL17_resp', 'IL23_resp', 'TNFa_resp',
        'IL17_non', 'IL23_non', 'TNFa_non'
    ])

    drug_targets = ['IL17', 'IL23', 'TNFa']

    # Combine all drug targets per endotype
    df['Total_resp'] = df[[f'{d}_resp' for d in drug_targets]].sum(axis=1)
    df['Total_non'] = df[[f'{d}_non' for d in drug_targets]].sum(axis=1)
    df['Total'] = df['Total_resp'] + df['Total_non']

    # Normalize totals
    df['Resp_prop'] = df['Total_resp'] / df['Total']
    df['Non_prop'] = df['Total_non'] / df['Total']

    # Compute p-value, OR and 95% CI vs all other endotypes
    stats_df = get_stats(df=df)
    stats_df.to_excel(os.path.join(save_folder, "Responder_distribution_stats.xlsx"))

    x = np.arange(len(df['Endotype']))

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    # -------- Plot: total responders vs non-responders --------
    ax.bar(x, df['Resp_prop'], label='Responders', color='green')
    ax.bar(x, df['Non_prop'], bottom=df['Resp_prop'], label='Non-responders', color='red', alpha=0.7)

    # # Add OR text above each bar
    # for i, row in stats_df.iterrows():
    #     ax.text(x[i], 1.02, f"p-value={row['p_value']}; OR={row['Odds_Ratio']:.2f}\n95%CI=[{row['CI_low']:.2f},{row['CI_high']:.2f}]",
    #             ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(df['Endotype'])
    ax.set_ylabel('Proportion of patients')
    ax.set_title('All drug targets combined')
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=False, shadow=False, ncol=5, frameon=False)

    sns.despine(fig=fig)

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "Responder_distribution.pdf"), bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Drug_target_independence', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)
    main(save_folder=output_dir)
