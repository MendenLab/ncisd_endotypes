from scripts.utils import add_colors
from scripts.feature_engineering import scaling
from scripts.utils import normalisation

import numpy as np
import scanpy as sc
import pandas as pd
from datetime import date
from scipy.spatial.distance import euclidean, pdist, squareform

import itertools

import os
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import pickle

fileformat = ".pdf"


def similarity_func(u, v):
    return 1 / (1 + euclidean(u, v))


def plot_umap_deltapga_scores(
    adata,
    observation_name,
    label,
    save_folder,
    use_rep="X_umap",
    color=None,
    figsize=(6, 6),
):
    mask_nonresponder = (adata.obs["Response week 12"] == "Non-responder") & (
        adata.obs["Therapeutic target"] == "IL-23"
    )
    mask_responder = (adata.obs["Response week 12"] == "Responder") & (
        adata.obs["Therapeutic target"] == "IL-23"
    )
    mask_superresponder = (adata.obs["Response week 12"] == "Super-responder") & (
        adata.obs["Therapeutic target"] == "IL-23"
    )

    fig, ax = plt.subplots(figsize=figsize)
    for ind_label, cluster_id in enumerate(label):
        inds = np.where(np.asarray(label) == cluster_id)[0]
        mask_cluster = adata.obs[observation_name] == cluster_id

        if len(inds) > 0:
            ax.scatter(
                adata.obsm[use_rep][:, 0][mask_cluster],
                adata.obsm[use_rep][:, 1][mask_cluster],
                marker=".",
                label=cluster_id,
                c=color[inds[0]],
                s=160,
            )
    # ax.legend(title=observation_name, bbox_to_anchor=(1, 1.02), loc='upper left', ncol=1, prop={'size': 10},
    #           frameon=False, title_fontsize=15)

    # non-responder
    ax.plot(
        adata.obsm[use_rep][:, 0][mask_nonresponder],
        adata.obsm[use_rep][:, 1][mask_nonresponder],
        marker="^",
        color="orangered",
        label="Non Responder",
        markersize=10,
        markeredgecolor="black",
        linestyle="none",
    )
    # responder
    ax.plot(
        adata.obsm[use_rep][:, 0][mask_responder],
        adata.obsm[use_rep][:, 1][mask_responder],
        marker="^",
        color="gold",
        label="Responder",
        markersize=10,
        markeredgecolor="black",
        linestyle="none",
    )
    # super responder
    ax.plot(
        adata.obsm[use_rep][:, 0][mask_superresponder],
        adata.obsm[use_rep][:, 1][mask_superresponder],
        marker="^",
        color="limegreen",
        label="Super Responder",
        markersize=10,
        markeredgecolor="black",
        linestyle="none",
    )

    ax.set_ylabel("UMAP2", fontsize=18)
    ax.set_xlabel("UMAP1", fontsize=18)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    plt.tight_layout()
    fig.savefig(
        os.path.join(
            save_folder,
            "UMAP_Response week 12_{}{}".format(observation_name, fileformat),
        )
    )
    plt.close(fig=fig)


def export_legend(legend, save_folder, ncol, key, filename="legend", expand=None):
    # https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib
    if expand is None:
        expand = [-5, -5, 5, 5]
    fig = legend.figure
    sns.despine(left=True, bottom=True, fig=fig)
    fig.canvas.draw()
    # bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    # fig.tight_layout()
    fig.savefig(
        os.path.join(save_folder, "{}_{}_{}{}".format(filename, ncol, key, fileformat)),
        dpi="figure",
        bbox_inches=bbox,
    )
    plt.close(fig=fig)


def load_molecular_subtypes(adata, data_root):
    try:
        with open(
            os.path.join(data_root, "Optimalres_0.5__OptimalGPTK_7.pkl"), "rb"
        ) as ff:
            df_res05 = pickle.load(ff)
        subtypes_res05 = df_res05.iloc[0]["Molecular subtypes"]
        with open(
            os.path.join(data_root, "Optimalres_0.6__OptimalGPTK_7.pkl"), "rb"
        ) as ff:
            df_res06 = pickle.load(ff)
        subtypes_res06 = df_res06.iloc[0]["Molecular subtypes"]
        with open(
            os.path.join(data_root, "Optimalres_0.7__OptimalGPTK_7.pkl"), "rb"
        ) as ff:
            df_res07 = pickle.load(ff)
        subtypes_res07 = df_res07.iloc[0]["Molecular subtypes"]
        with open(
            os.path.join(data_root, "Optimalres_0.8__OptimalGPTK_7.pkl"), "rb"
        ) as ff:
            df_res08 = pickle.load(ff)
        subtypes_res08 = df_res08.iloc[0]["Molecular subtypes"]
        with open(
            os.path.join(data_root, "Optimalres_0.9__OptimalGPTK_7.pkl"), "rb"
        ) as ff:
            df_res09 = pickle.load(ff)
        subtypes_res09 = df_res09.iloc[0]["Molecular subtypes"]
    except ModuleNotFoundError:
        df_res05 = pd.read_pickle(
            os.path.join(data_root, "Optimalres_0.5__OptimalGPTK_7.pkl")
        )
        subtypes_res05 = df_res05.iloc[0]["Molecular subtypes"]
        df_res06 = pd.read_pickle(
            os.path.join(data_root, "Optimalres_0.6__OptimalGPTK_7.pkl")
        )
        subtypes_res06 = df_res06.iloc[0]["Molecular subtypes"]
        df_res07 = pd.read_pickle(
            os.path.join(data_root, "Optimalres_0.7__OptimalGPTK_7.pkl")
        )
        subtypes_res07 = df_res07.iloc[0]["Molecular subtypes"]
        df_res08 = pd.read_pickle(
            os.path.join(data_root, "Optimalres_0.8__OptimalGPTK_7.pkl")
        )
        subtypes_res08 = df_res08.iloc[0]["Molecular subtypes"]
        df_res09 = pd.read_pickle(
            os.path.join(data_root, "Optimalres_0.9__OptimalGPTK_7.pkl")
        )
        subtypes_res09 = df_res09.iloc[0]["Molecular subtypes"]

    assert np.all(
        df_res05.iloc[0]["MUC IDs"] == adata.obs.index
    ), "Index are not in the same order as in adata"

    adata.obs["Molecular Subtype res0.5"] = subtypes_res05 + 1
    adata.obs["Molecular Subtype res0.5"] = adata.obs[
        "Molecular Subtype res0.5"
    ].astype("category")
    adata.obs["Molecular Subtype res0.6"] = subtypes_res06 + 1
    adata.obs["Molecular Subtype res0.6"] = adata.obs[
        "Molecular Subtype res0.6"
    ].astype("category")
    adata.obs["Molecular Subtype res0.7"] = subtypes_res07 + 1
    adata.obs["Molecular Subtype res0.7"] = adata.obs[
        "Molecular Subtype res0.7"
    ].astype("category")
    adata.obs["Molecular Subtype res0.8"] = subtypes_res08 + 1
    adata.obs["Molecular Subtype res0.8"] = adata.obs[
        "Molecular Subtype res0.8"
    ].astype("category")
    adata.obs["Molecular Subtype res0.9"] = subtypes_res09 + 1
    adata.obs["Molecular Subtype res0.9"] = adata.obs[
        "Molecular Subtype res0.9"
    ].astype("category")

    # Add embedding
    with open(os.path.join(data_root, "MolecularSubtypes_embedding.pkl"), "rb") as ff:
        embedding_geneselection = pickle.load(ff)

    adata.obsm["X_geneselection_umap"] = np.asarray(embedding_geneselection)
    adata.var["geneselection"] = adata.var_names.isin(df_res09.iloc[0]["Gene names"])

    return adata


def plot_barplot_counts(adata, save_folder, obs, title):
    df = adata.obs[[obs, "Endotypes"]]
    df_plot = (
        df.groupby([obs, "Endotypes"])
        .size()
        .reset_index()
        .pivot(columns=obs, index="Endotypes", values=0)
    )

    fig, ax = plt.subplots()
    df_plot.plot.bar(
        rot=0,
        stacked=True,
        color=adata.uns["{}_colors".format(obs)],
        figsize=(8, 6),
        width=0.95,
        edgecolor="black",
        linewidth=0.1,
        ax=ax,
        legend=False,
    )
    ax.set_xlabel("")
    ax.set_ylabel("counts", fontsize=18)
    ax.spines["left"].set_visible(True)
    ax.tick_params(axis="both", labelsize=15)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(True)

    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        if height >= 0.05:
            ax.text(
                x + width / 2,
                y + height / 2,
                "{:.0f}".format(height),
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=12,
            )

    fig.legend(
        title=title,  # obs.capitalize()
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        ncol=3,
        prop={"size": 15},
        frameon=False,
        title_fontsize=17,
    )

    plt.savefig(
        os.path.join(save_folder, "Barplot_Endotypes_vs_{}_counts.pdf".format(obs)),
        bbox_inches="tight",
    )
    plt.close()


def main(save_folder):
    data_root = os.path.join(
        "/Volumes/CH__data/Projects/Eyerich_AG_projects",
        "BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/input/h5_files/LESION",
    )
    adata = sc.read(
        os.path.join(
            data_root,
            "__".join(["Lesion_RNAseq_20210720_patient_meta_data_v04",
                       "CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected", "Endotypes_230620.h5"]),
        )
    )
    # adata = load_molecular_subtypes(adata=adata, data_root=data_root)

    for obs_name in [
        "Molecular Subtype res0.5",
        "Molecular Subtype res0.6",
        "Molecular Subtype res0.7",
        "Molecular Subtype res0.8",
        "Molecular Subtype res0.9",
        'Endotypes'
    ]:
        if obs_name != "Molecular Subtype res0.9":
            adata, _ = add_colors.endotype_order_lesion(adata=adata, obs_name=obs_name)

        ax = sns.countplot(
            x=obs_name, data=adata.obs, palette=adata.uns["{}_colors".format(obs_name)]
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.savefig(os.path.join(save_folder, "Countplot_{}.pdf".format(obs_name)))
        plt.close("all")

    # Number of columns in legend
    ncol = 1
    use_rep = "X_geneselection_umap"

    pd.crosstab(
        adata.obs["Pattern"], adata.obs["Molecular Subtype res0.5"], dropna=False
    )
    pd.crosstab(adata.obs["diag"], adata.obs["Molecular Subtype res0.5"], dropna=False)

    plot_barplot_counts(
        adata=adata, save_folder=save_folder, obs="diag", title="diagnosis"
    )
    plot_barplot_counts(
        adata=adata, save_folder=save_folder, obs="Pattern", title="pattern"
    )

    adata.obs["MUC IDs"] = adata.obs.index
    adata.obs.to_csv(os.path.join(save_folder, "Metadata.csv"))

    # Get normed and scaled counts
    bulk_data_temp = adata.copy()
    df_gex = bulk_data_temp.to_df(layer="counts")
    # Apply normalisation and batch correction
    array_gex = normalisation.edger_normalise(
        X=np.asarray(df_gex), coldata=bulk_data_temp.obs, batch_keys=['batchID', 'sex'])
    df_gex_normed = pd.DataFrame(array_gex, columns=df_gex.columns, index=df_gex.index)
    # Scale data on selected genes
    df_gex_scaled, _ = scaling.scale_gex(data_train=df_gex_normed, sparse_data=False)

    genes_of_interest = [
        "IL6",
        "TNF",
        "IL26",
        "IL23A",
        "IL17A",
        "MS4A1",
        "IL31",
        "IL13",
        "IL4R",
        "IL17F",
    ]
    for gene in genes_of_interest:
        # counts = adata.to_df(layer='counts')[gene]
        counts = df_gex_scaled[gene]

        vcenter = 0
        vmin, vmax = counts.min(), counts.max()
        normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
        colormap = sns.color_palette("coolwarm", as_cmap=True)

        fig, ax = plt.subplots()
        sns.scatterplot(
            y=adata.obsm[use_rep][:, 1],
            x=adata.obsm[use_rep][:, 0],
            c=counts,
            linewidth=0,
            norm=normalize,
            cmap=colormap,
            ax=ax,
        )
        scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
        scalarmappaple.set_array(counts)
        cbar = fig.colorbar(scalarmappaple, ax=ax)
        cbar.ax.get_yaxis().labelpad = 10
        cbar.ax.set_ylabel("scaled counts", rotation=90, fontsize=14)
        ax.set_title(gene, fontsize=18)
        ax.set_ylabel("UMAP2", fontsize=18)
        ax.set_xlabel("UMAP1", fontsize=18)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        plt.savefig(
            os.path.join(
                save_folder, "Scanpy_UMAP_{}_userep{}.pdf".format(gene, use_rep)
            )
        )
        plt.close(fig=fig)

        # fig, ax = plt.subplots()
        # im = ax.scatter(adata.obsm[use_rep][:, 0], adata.obsm[use_rep][:, 1],  marker='.', c=counts, s=160,
        #                 cmap=sns.color_palette("coolwarm", as_cmap=True))
        # cbar = fig.colorbar(im, ax=ax)
        # cbar.ax.get_yaxis().labelpad = 10
        # cbar.ax.set_ylabel('scaled counts', rotation=90, fontsize=14)
        # ax.set_title(gene, fontsize=18)
        # ax.set_ylabel('UMAP2', fontsize=18)
        # ax.set_xlabel('UMAP1', fontsize=18)
        # ax.spines["top"].set_visible(False)
        # ax.spines["right"].set_visible(False)
        # fig.tight_layout()
        # plt.savefig(os.path.join(save_folder, "Scanpy_UMAP_{}_userep{}.pdf".format(gene, use_rep)))
        # plt.close(fig=fig)

    adata.obs["empty"] = "unknown"
    adata.obs["empty"] = adata.obs["empty"].astype("category")
    adata.uns["empty_colors"] = ["grey"]

    adata.uns["{}_colors".format('Diag_PGA')] = ['green', 'blue', 'orange', 'red']

    for obs in [
        "diag",
        "sdiag",
        "Pattern",
        "Molecular Subtype res0.5",
        "Molecular Subtype res0.6",
        "Molecular Subtype res0.7",
        "Molecular Subtype res0.8",
        "Molecular Subtype res0.9",
        "Endotypes",
        "empty",
        'Diag_PGA'
    ]:
        labels = list(adata.obs[obs].cat.categories)
        color = list(adata.uns["{}_colors".format(obs)])
        fig, ax = plt.subplots(figsize=(6, 6))
        for ind_label, cluster_id in enumerate(labels):
            inds = np.where(np.asarray(labels) == cluster_id)[0]
            mask_cluster = adata.obs[obs] == cluster_id

            if len(inds) > 0:
                ax.scatter(
                    adata.obsm[use_rep][:, 0][mask_cluster],
                    adata.obsm[use_rep][:, 1][mask_cluster],
                    marker=".",
                    label=cluster_id,
                    c=color[inds[0]],
                    s=160,
                )
        ax.set_ylabel("UMAP2", fontsize=18)
        ax.set_xlabel("UMAP1", fontsize=18)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        plt.savefig(
            os.path.join(
                save_folder, "Scanpy_UMAP_{}_userep{}.pdf".format(obs, use_rep)
            )
        )
        plt.close(fig=fig)

        plot_umap_deltapga_scores(
            adata=adata,
            observation_name=obs,
            label=list(adata.obs[obs].cat.categories),
            save_folder=save_folder,
            use_rep=use_rep,
            color=adata.uns["{}_colors".format(obs)],
            figsize=(6, 6),
        )

        # Plot legend separately
        def f(c): return plt.plot([], [], marker="o", color=c, ls="none")[0]
        handles = [f(c=color[i]) for i in range(len(color))]

        legend = plt.legend(
            handles, labels, loc=3, framealpha=1, frameon=False, ncol=ncol
        )
        export_legend(
            legend, save_folder=save_folder, ncol=ncol, key=obs, expand=[-5, -5, 5, 5]
        )

    use_rep = "X_umap"
    for obs in [
        "diag",
        "sdiag",
        "Pattern",
        "Molecular Subtype res0.5",
        "Molecular Subtype res0.6",
        "Molecular Subtype res0.7",
        "Molecular Subtype res0.8",
        "Molecular Subtype res0.9",
        "Endotypes",
        "empty",
        'Diag_PGA'
    ]:
        labels = list(adata.obs[obs].cat.categories)
        if "{}_colors".format(obs) in adata.uns.keys():
            color = list(adata.uns["{}_colors".format(obs)])
        else:
            color = sc.pl.palettes.default_28[: len(labels)]
            adata.uns["{}_colors".format(obs)] = color
        fig, ax = plt.subplots(figsize=(6, 6))
        for ind_label, cluster_id in enumerate(labels):
            inds = np.where(np.asarray(labels) == cluster_id)[0]
            mask_cluster = adata.obs[obs] == cluster_id

            if len(inds) > 0:
                ax.scatter(
                    adata.obsm[use_rep][:, 0][mask_cluster],
                    adata.obsm[use_rep][:, 1][mask_cluster],
                    marker=".",
                    label=cluster_id,
                    c=color[inds[0]],
                    s=160,
                )
        ax.set_ylabel("UMAP2", fontsize=18)
        ax.set_xlabel("UMAP1", fontsize=18)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        plt.savefig(
            os.path.join(
                save_folder, "Scanpy_UMAP_{}_userep{}.pdf".format(obs, use_rep)
            ),  bbox_inches='tight'
        )
        plt.close(fig=fig)

    # Calculate similarity between clusters by average distance between every point and every group of points
    # 1. Calculate distance matrix
    dist_ma = pdist(df_gex_normed, similarity_func)
    # Convert to square form
    df_euclid = pd.DataFrame(
        squareform(dist_ma), columns=df_gex_normed.index, index=df_gex_normed.index
    )
    df_euclid["cluster"] = adata.obs["Molecular Subtype res0.9"]
    # Read out per pairwise clusters the similarity values -> get average?
    pair_order_list = list(
        itertools.permutations(
            list(adata.obs["Molecular Subtype res0.9"].cat.categories), 2
        )
    )
    sim = []
    for pair in pair_order_list:
        names = list(
            df_euclid.loc[df_euclid["cluster"].isin([pair[0], pair[1]]), :].index
        )
        # df_euclid.loc[df_euclid['cluster'].isin(pair[0], pair[1]), names]
        sim_tmp = np.mean(
            df_euclid.loc[df_euclid["cluster"].isin([pair[0], pair[1]]), names]
        )
        sim.append(sim_tmp)

        print("SIM {} vs {}: {}".format(pair[0], pair[1], sim_tmp))

    # TODO Create boxplots showing with mean expression of all genes per sample for each cluster
    df_gex_normed_mean = df_gex_normed.loc[:, df_gex_normed.columns != "cluster"].mean(
        axis=1
    )
    df_gex_normed_mean["cluster"] = adata.obs["Molecular Subtype res0.9"]
    df_gex_normed_mean.groupby("cluster").boxplot()

    print("Done")


if __name__ == "__main__":
    output_dir = os.path.join(
        "/Volumes",
        "CH__data",
        "Projects",
        "Eyerich_AG_projects",
        "BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer",
        "analysis",
        "Molecular_subtypes",
        "output",
        "Figure_2I",
        str(date.today()),
    )
    os.makedirs(output_dir, exist_ok=True)
    main(save_folder=output_dir)
