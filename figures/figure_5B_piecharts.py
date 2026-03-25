import pandas as pd
from datetime import date
import os

import matplotlib.pyplot as plt


def main(save_folder):
    excel_path = "/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/raw_data/Table_S9.xlsx"  # path to your Excel file
    endotypes_to_plot = ["E5", "E11", "E12", "E13"]

    # Category labels and colors
    category_labels = {
        "SR": "Super Responder",
        "R": "Responder",
        "NR": "Non-Responder",
        "SNR": "Super Non-Responder"
    }
    colors = {
        "Super Responder": "#006400",      # dark green
        "Responder": "#90EE90",            # light green
        "Non-Responder": "#FF9999",        # light red
        "Super Non-Responder": "#8B0000"   # dark red
    }
    # Define order for pie slices
    category_order = ["SNR", "NR", "R", "SR"]

    # === READ ALL SHEETS ===
    sheets = pd.read_excel(excel_path, sheet_name=None, index_col=0)

    # Skip first sheet
    sheet_items = list(sheets.items())[1:]

    for sheet_name, df in sheet_items:
        print(f"Processing sheet: {sheet_name}")

        # Split drug and category
        df = df.copy()
        df["Drug"] = df.index.str.split().str[0]
        df["Category"] = df.index.str.split().str[1]

        for drug, drug_df in df.groupby("Drug"):
            fig, axes = plt.subplots(1, len(endotypes_to_plot), figsize=(14, 4))

            # Ensure all categories exist (fill missing with 0)
            drug_df_full = pd.DataFrame({"Category": category_order})
            drug_df_full["Drug"] = drug
            drug_df_full = pd.merge(
                drug_df_full,
                drug_df[["Category"] + endotypes_to_plot],
                on="Category",
                how="left"
            ).fillna(0)

            for ax, endotype in zip(axes, endotypes_to_plot):
                counts = drug_df_full[endotype].values
                labels = [category_labels[cat] for cat in drug_df_full["Category"]]
                cols = [colors[label] for label in labels]

                # Remove categories with 0 counts to keep pie clean
                mask = counts > 0
                counts = counts[mask]
                labels = [l for l, m in zip(labels, mask) if m]
                cols = [c for c, m in zip(cols, mask) if m]

                ax.pie(counts, labels=labels, autopct="%1.0f%%", colors=cols, startangle=90)
                ax.set_title(f"{drug} - {endotype}")

            fig.suptitle(f"{sheet_name} - {drug} Response Distribution", fontsize=14)
            plt.tight_layout()
            plt.savefig(os.path.join(save_folder, '{}_Piecharts_{}.pdf'.format(sheet_name, drug)),
                        bbox_inches='tight')
            plt.close(fig=fig)


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Drug_responds_piecharts', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
