import pandas as pd
import numpy as np
import math
import re
import os
from datetime import date

from scipy.stats import fisher_exact, norm

import matplotlib.pyplot as plt
import seaborn as sns


def odds_ratio_ci(a, b, c, d, alpha=0.05):
    # Add 0.5 continuity correction (Haldane-Anscombe correction) if any cell is zero
    if 0 in [a, b, c, d]:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
    or_val = (a * d) / (b * c)
    log_or = math.log(or_val)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    z = norm.ppf(1 - alpha / 2)
    lower = math.exp(log_or - z * se)
    upper = math.exp(log_or + z * se)
    return or_val, lower, upper


def cohens_h(p1, p2):
    p1, p2 = np.clip(p1, 1e-9, 1-1e-9), np.clip(p2, 1e-9, 1-1e-9)
    return 2 * (np.arcsin(np.sqrt(p1)) - np.arcsin(np.sqrt(p2)))


def get_annotations(df, endotypes, drug_targets):
    # For collecting annotation lines per endotype
    annotations_dict = {e: [] for e in endotypes}

    # Calculate OR, CI, p, Cohen's h for each endotype and drug target pair
    for idx, row in df.iterrows():
        endotype = row['Endotype']
        ann_lines = []
        for i in range(len(drug_targets)):
            for j in range(i+1, len(drug_targets)):
                dt1, dt2 = drug_targets[i], drug_targets[j]
                a = row[f'{dt1}_resp']
                b = row[f'{dt1}_non']
                c = row[f'{dt2}_resp']
                d = row[f'{dt2}_non']
                table = [[a, b], [c, d]]
                _, p_value = fisher_exact(table, alternative='greater')
                or_val, l_ci, u_ci = odds_ratio_ci(a, b, c, d)

                p1 = a / (a + b)
                p2 = c / (c + d)
                h = cohens_h(p1, p2)

                # Short annotation
                ann_lines.append(f"{dt1} vs {dt2}: p={p_value:.2f} | h={h:.2f} | OR={or_val:.2f} [{l_ci:.2f}-{u_ci:.2f}]")
        annotations_dict[endotype] = ann_lines

    return annotations_dict


def get_stats_dataframe(df, annotations_dict):
    rows = []

    for endotype, annotations in annotations_dict.items():
        # get the row from df corresponding to this endotype
        row = df[df['Endotype'] == endotype].iloc[0]
        for ann in annotations:
            comp_match = re.match(r'(.+):', ann)
            comparison = comp_match.group(1) if comp_match else ''

            # extract the two drug targets (e.g., "IL17 vs IL23")
            dt1, dt2 = [s.strip() for s in comparison.split('vs')]

            p_match = re.search(r'p=([\d\.]+)', ann)
            p_value = float(p_match.group(1)) if p_match else None

            h_match = re.search(r'h=([\d\.\-]+)', ann)
            cohens_h = float(h_match.group(1)) if h_match else None

            or_match = re.search(r'OR=([\d\.]+) \[([\d\.]+)-([\d\.]+)\]', ann)
            if or_match:
                odds_ratio = float(or_match.group(1))
                ci_lower = float(or_match.group(2))
                ci_upper = float(or_match.group(3))
            else:
                odds_ratio = ci_lower = ci_upper = None

            # Get responder/non-responder counts from original df row
            resp_dt1 = row[f'{dt1}_resp']
            nonresp_dt1 = row[f'{dt1}_non']
            resp_dt2 = row[f'{dt2}_resp']
            nonresp_dt2 = row[f'{dt2}_non']

            rows.append({
                'Endotype': endotype,
                'Comparison': comparison,
                'p_value': p_value,
                'Cohens_h': cohens_h,
                'Odds_Ratio': odds_ratio,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper,
                f'{dt1}_resp': resp_dt1,
                f'{dt1}_nonresp': nonresp_dt1,
                f'{dt2}_resp': resp_dt2,
                f'{dt2}_nonresp': nonresp_dt2
            })

    annotations_df = pd.DataFrame(rows)
    return annotations_df


def main(df, save_folder):
    drug_targets = ['IL17', 'IL23', 'TNFa']
    colors = {'IL17': '#1f77b4', 'IL23': '#ff7f0e', 'TNFa': '#2ca02c'}
    # colors = {'IL17': 'red', 'IL23': 'blue', 'TNFa': 'purple'}

    endotypes = df['Endotype'].tolist()
    bar_width = 0.25
    x = np.arange(len(endotypes))

    annotations_dict = get_annotations(df=df, endotypes=endotypes, drug_targets=drug_targets)
    # Create DataFrame with all annotations
    annotations_df = get_stats_dataframe(df=df, annotations_dict=annotations_dict)
    annotations_df.to_excel(os.path.join(save_folder, 'Drugtarget_ranking_stats.xlsx'))

    fig, ax = plt.subplots(figsize=(8, 4))
    for i, dt in enumerate(drug_targets):
        prop_resp = []
        prop_nonresp = []
        for _, row in df.iterrows():
            resp = row[f'{dt}_resp']
            nonresp = row[f'{dt}_non']
            total = resp + nonresp
            prop_resp.append(resp / total)
            prop_nonresp.append(nonresp / total)

        prop_resp = np.array(prop_resp)
        prop_nonresp = np.array(prop_nonresp)

        bars1 = ax.bar(x + i * bar_width, prop_resp, width=bar_width, color=colors[dt], label=f'{dt} responders')
        bars2 = ax.bar(x + i * bar_width, prop_nonresp, bottom=prop_resp, width=bar_width, color=colors[dt], alpha=0.5,
                       label=f'{dt} non-responders')

        # Add percentage labels on each segment (responders)
        for rect, val in zip(bars1, prop_resp):
            if val > 0:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width() / 2, height / 2 + rect.get_y(),
                        f'{val*100:.1f}%', ha='center', va='center', color='black', fontsize=6, fontweight='bold')

        # Add percentage labels on each segment (non-responders)
        for rect, val, bottom in zip(bars2, prop_nonresp, prop_resp):
            if val > 0:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width() / 2, bottom + height / 2,
                        f'{val*100:.1f}%', ha='center', va='center', color='black', fontsize=6, fontweight='bold')

    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(endotypes)
    ax.set_ylabel('Proportion')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left')

    ax.axhline(0.5, color='gray', linestyle='--', linewidth=1)

    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()
    sns.despine(fig=fig)
    plt.savefig(os.path.join(save_folder, "Responder_enrichment.pdf"), bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    key = 'SR'  #'SR'  #'SR+R'

    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Drug_target_ranking_{}'.format(key), str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    # if key == 'SR+R':
    #     # combines responders (SR + R) and non-responders (NR + SNR)
    #     df_counts_infos = pd.DataFrame([
    #         ['E5',  ['E5', 'E11/E12/E13'],  5, 5, 6, 0, 1, 4],
    #         ['E11', ['E11', 'E5/E12/E13'],  6, 7, 11, 1, 7, 14],
    #         ['E12', ['E12', 'E5/E11/E13'],  8, 27, 15, 3, 7, 17],
    #         ['E13', ['E13', 'E5/E11/E12'],  5, 4, 8, 1, 2, 8],
    #     ], columns=[
    #         'Endotype', 'Others',
    #         'IL17_resp', 'IL23_resp', 'TNFa_resp',
    #         'IL17_non', 'IL23_non', 'TNFa_non'
    #     ])
    # else:
    #     # considers only super-responders (SR) and non-responders (NR + SNR)
    #     df_counts_infos = pd.DataFrame([
    #         ['E5',  ['E5', 'E11/E12/E13'],  2, 2, 0, 0, 1, 4],  #  'IL17_resp', 'IL23_resp', 'TNFa_resp', 'IL17_non', 'IL23_non', 'TNFa_non'
    #         ['E11', ['E11', 'E5/E12/E13'],  4, 2, 3, 1, 7, 14],
    #         ['E12', ['E12', 'E5/E11/E13'],  4, 5, 0, 3, 7, 17],
    #         ['E13', ['E13', 'E5/E11/E12'],  2, 1, 2, 1, 2, 8],
    #     ], columns=[
    #         'Endotype', 'Others',
    #         'IL17_resp', 'IL23_resp', 'TNFa_resp',
    #         'IL17_non', 'IL23_non', 'TNFa_non'
    #     ])

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
