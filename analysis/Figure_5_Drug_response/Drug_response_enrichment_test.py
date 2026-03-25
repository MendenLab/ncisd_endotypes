import numpy as np
from scipy.stats import fisher_exact
import math
from scipy.stats import norm
import pandas as pd
from datetime import date
import os

import matplotlib.pyplot as plt
import seaborn as sns


def plot_barchart(responders, non_responders, group_name, p_value, cohens_h_value, odds_ratio, l_ci, u_ci,
                  save_folder, endotype, drug_target):
    # Plot: Normalized (100%) Stacked Bar Chart -> Visualised proportions, tied to statistcal testing
    # Convert to proportions
    totals = [r + n for r, n in zip(responders, non_responders)]
    prop_responders = [r / t for r, t in zip(responders, totals)]
    prop_non_responders = [nr / t for nr, t in zip(non_responders, totals)]

    fig, ax = plt.subplots(figsize=(6, 4))

    bars1 = ax.bar(group_name, prop_responders, label='Responders')
    bars2 = ax.bar(group_name, prop_non_responders, bottom=prop_responders, label='Non-responders')

    # Annotate with proportions
    for i in range(len(group_name)):
        ax.text(i, prop_responders[i]/2, f"{prop_responders[i]*100:.1f}%",
                ha='center', va='center', color='white', fontsize=10)
        ax.text(i, prop_responders[i] + prop_non_responders[i]/2,
                f"{prop_non_responders[i]*100:.1f}%", ha='center', va='center', color='white', fontsize=10)

    # Optional annotation for p-value or effect size
    ax.text(0.5, 1.05,
            f"p = {p_value:.2f} | Cohen's h = {cohens_h_value:.2f} | OR = {odds_ratio:.2f} (95% CI: {l_ci:.2f}–{u_ci:.2f})",
            ha='center', fontsize=10)

    # Labels and formatting
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Proportion")
    ax.set_title("Proportion of Responders by Group (Normalized)")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
              ncol=3, fancybox=False, shadow=False, frameon=False)
    sns.despine(fig=fig)

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "{}_{}_responder_enrichment.pdf".format(endotype, drug_target)),
                bbox_inches='tight')
    plt.close()


# Cohen's h formula
def cohens_h(p1, p2):
    return 2 * (math.asin(math.sqrt(p1)) - math.asin(math.sqrt(p2)))


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


def main(df, save_folder):
    df_stats = pd.DataFrame()
    for endotype in df['Endotype']:
        group_name = df.loc[df['Endotype'] == endotype, 'Others'].item()
        for drug_target in ['IL23', 'IL17', 'TNFa']:
            responder_name = '{}_resp'.format(drug_target)
            nonresponder_name = '{}_non'.format(drug_target)

            responders_target = df.loc[df['Endotype'] == endotype, responder_name].item()
            total_events_target = responders_target + df.loc[df['Endotype'] == endotype, nonresponder_name].item() # includes non-responders

            responders_others = df.loc[df['Endotype'] != endotype, responder_name].sum()
            total_events_others = responders_others + df.loc[df['Endotype'] != endotype, nonresponder_name].sum()

            non_responders = [total_events_target - responders_target, total_events_others - responders_others]

            # 2x2 contingency table
            #           Responders  Non-responders
            # table = [[27, 7],       # E12
            #          [16, 10]]       # E5/E11/E13
            table = [[responders_target, non_responders[0]],
                     [responders_others, non_responders[1]]]

            # Fisher's Exact Test (default alternative='two-sided')
            # Add 0.5 continuity correction (modified Haldane-Anscombe correction) if any cell is zero
            if 0 in [table[0][0], table[0][1], table[1][0], table[1][1]]:
                table[0][0] += 0.5
                table[0][1] += 0.5
                table[1][0] += 0.5
                table[1][1] += 0.5
            odds_ratio, p_value = fisher_exact(table, alternative='greater')

            print("Fisher's Exact Test ({} > {}): p = {:.4f}".format(group_name[0], group_name[1], p_value))  # E12: 0.11
            print(f"Odds Ratio = {odds_ratio:.2f}")  # Odds Ratio = 2.41 → more likely to respond in E12
            # Odds Ratio (OR):
            # OR > 1 → more likely to respond in E12
            # OR = 1 → no difference
            # OR < 1 → less likely to respond in E12
            # Calculate CI of OR
            l_ci, u_ci = odds_ratio_ci(table=table, odds_ratio=odds_ratio, alpha=0.05)

            # Proportions
            # Within-group proportions
            p1_within = responders_target / total_events_target
            p2_within = responders_others / total_events_others

            # Overall proportions (relative to total N)
            # p1_total = responders_target / total_all_patients
            # p2_total = responders_others / total_all_patients

            h = cohens_h(p1_within, p2_within)

            print(f"Proportion {group_name[0]}: {p1_within:.3f}")
            print(f"Proportion {group_name[1]}: {p2_within:.3f}")
            print(f"Cohen's h = {h:.3f}")  # E12 Cohen's h = 0.396 → small effect
            # h ≈ 0.2 → small effect
            # h ≈ 0.5 → medium effect
            # h ≈ 0.8+ → large effect

            # Plot
            plot_barchart(responders=[responders_target, responders_others], non_responders=non_responders,
                          group_name=group_name, p_value=p_value, l_ci=l_ci, u_ci=u_ci,
                          cohens_h_value=h, odds_ratio=odds_ratio, save_folder=save_folder, endotype=endotype,
                          drug_target=drug_target)

            # Store in dataframe
            df_tmp = pd.DataFrame([{
                'Drug_target': drug_target,
                "Group1": group_name[0],
                "Group2": group_name[1],
                "Responders_Group1": responders_target,
                "NonResponders_Group1": non_responders[0],
                "Responders_Group2": responders_others,
                "NonResponders_Group2": non_responders[1],
                "OddsRatio": odds_ratio,
                "CI_lower": l_ci,
                "CI_upper": u_ci,
                "p_value": p_value,
                "Prop_within_Group1": p1_within,
                "Prop_within_Group2": p2_within,
                "Cohen_h": h
            }])

            df_stats = pd.concat([df_stats, df_tmp])
    df_stats.to_excel(os.path.join(save_folder, 'Drugresponse_tests.xlsx'), index=False)


if __name__ == '__main__':
    key = 'SR'  #'SR'  #'SR+R'
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Drug_responds_test_{}'.format(key), str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    if key == 'SR+R':
        # combines responders (SR + R) and non-responders (NR + SNR)
        df_counts_infos = pd.DataFrame([
            ['E5',  ['E5', 'E11/E12/E13'],  3, 3, 6, 0, 1, 4], #  'IL17_resp', 'IL23_resp', 'TNFa_resp', 'IL17_non', 'IL23_non', 'TNFa_non'
            ['E11', ['E11', 'E5/E12/E13'],  5, 6, 9, 2, 7, 13],
            ['E12', ['E12', 'E5/E11/E13'],  5, 24, 15, 3, 7, 15],
            ['E13', ['E13', 'E5/E11/E12'],  4, 4, 7, 1, 2, 8],
        ], columns=[
            'Endotype', 'Others',
            'IL17_resp', 'IL23_resp', 'TNFa_resp',
            'IL17_non', 'IL23_non', 'TNFa_non'
        ])
    else:
        # considers only super-responders (SR) and non-responders (NR + SNR)
        df_counts_infos = pd.DataFrame([
            ['E5',  ['E5', 'E11/E12/E13'],  2, 2, 0, 0, 1, 4],  #  'IL17_resp', 'IL23_resp', 'TNFa_resp', 'IL17_non', 'IL23_non', 'TNFa_non'
            ['E11', ['E11', 'E5/E12/E13'],  4, 2, 3, 2, 7, 13],
            ['E12', ['E12', 'E5/E11/E13'],  4, 7, 0, 3, 7, 15],
            ['E13', ['E13', 'E5/E11/E12'],  2, 1, 2, 1, 2, 8],
        ], columns=[
            'Endotype', 'Others',
            'IL17_resp', 'IL23_resp', 'TNFa_resp',
            'IL17_non', 'IL23_non', 'TNFa_non'
        ])

    main(df=df_counts_infos, save_folder=output_dir)
