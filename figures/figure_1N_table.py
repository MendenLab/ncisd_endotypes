import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import numpy as np
from datetime import date
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

import os

matplotlib.rcParams['axes.axisbelow'] = True


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def gradient_fill(x, y, fill_color=None, ax=None, **kwargs):
    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.
    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.
    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.
    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """
    if ax is None:
        ax = plt.gca()

    line, = ax.plot(x, y, lw=0.0, **kwargs)
    if fill_color is None:
        fill_color = line.get_color()

    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 1.0 if alpha is None else alpha

    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:,:,:3] = rgb
    z[:,:,-1] = np.linspace(0, alpha, 100)[:,None]

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)

    xy = np.column_stack([x, y])
    xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    clip_path = Polygon(xy, facecolor='none', edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)

    ax.autoscale(True)
    return line, im


def gradient_image(ax, extent, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    extent
        The extent of the image as (xmin, xmax, ymin, ymax).
        By default, this is in Axes coordinates but may be
        changed using the *transform* kwarg.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular useful is *cmap*.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, extent=extent, interpolation='bicubic',
                   vmin=0, vmax=1, **kwargs)
    return im


def main(input_dir, save_folder):
    df = pd.read_excel(os.path.join(input_dir, 'Kilians_simplified_table_PS_v6.xlsx'),
                       sheet_name='Table_for_plot')
    df['efficacy'] = df['efficacy'].astype('category')
    df['value'] = df['value'].astype('category')
    df['color_efficacy'] = df['color_efficacy'].astype('category')
    df['disease'] = df['disease'].astype('category')
    df['drug'] = df['drug'].astype('category')

    # # difference Table v5 to v6
    # for col in df.columns:
    #     mask = df_v5[col] == df[col]
    #     if len(df.loc[~mask, col][pd.notna(df.loc[~mask, col])]) > 0:
    #         print('{}: '.format(col), df.loc[~mask, col][pd.notna(df.loc[~mask, col])])

    # df['value'].unique()
    # Add markers for efficacy level
    # markers = {'conflict': 'v', 'ineffective': 'o',
    #            'low': '^', 'moderate': '^', 'high': '^', 'very high': '^'}
    df['efficacy_markers'] = df['efficacy'].replace({
        'conflict': 'conflict', 'ineffective': 'ineffective', 'low': 'effective', 'moderate': 'effective',
        'high': 'effective', 'very high': 'effective'})
    df['efficacy_markers'] = df['efficacy_markers'].astype('category')
    df['efficacy_markers'] = df['efficacy_markers'].cat.reorder_categories(['effective', 'ineffective', 'conflict'])
    # triangle downwards ="ineffective" and circle for "conflict", triangle pointing upwards = "effective"
    markers = {'conflict': 'o', 'ineffective': 'v', 'effective': '^'}

    # Get colors for efficay level
    color_efficacy = {'conflict': 'darkblue', 'ineffective': 'lightgrey',
                      'low': 'lightgreen', 'moderate': 'mediumseagreen', 'high': 'forestgreen', 'very high': 'darkgreen'}

    df['marker_clinical_trial_ongoing'] = df['clinical_trial_ongoing']
    df.loc[df['marker_clinical_trial_ongoing'] == 'yes', 'marker_clinical_trial_ongoing'] = '*'
    df['marker_possible_induction'] = df['possible_induction']
    df.loc[df['marker_possible_induction'] == 'yes', 'marker_possible_induction'] = '-'
    df["efficacy_combined"] = df["marker_clinical_trial_ongoing"].astype(str) + df["marker_possible_induction"].astype(str)
    df['efficacy_combined'] = df['efficacy_combined'].astype(str).str.replace('nan', '')

    # Plot
    df['row_coord'] = df['target'].replace({'CD20': 4, 'IgE ': 3, 'IL-1': 9, 'IL-13': 1, 'IL-17': 5, 'IL-23': 6,
                                            'IL-31': 2, 'IL-36': 7, 'IL-4RA': 0, 'IL-6': 10, 'TNF-α': 8})
    # sort target:  IL-4RA, IL-13, IL-31, IgE , CD20, IL-17, IL-23, IL-36, TNF-α, IL-1, IL-6
    df['target'] = df['target'].astype('category')
    df['target'] = df['target'].cat.reorder_categories([
        'IL-4RA', 'IL-13', 'IL-31', 'IgE ', 'CD20', 'IL-17', 'IL-23', 'IL-36', 'TNF-α', 'IL-1', 'IL-6'])
    df['drug'] = df['drug'].cat.reorder_categories([
        'Dupilumab', 'Tralokinumab, Lebrikizumab', 'Nemolizumab', 'Omalizumab', 'Rituximab',
        'Bimekizumab, Brodalumab, Ixekizumab, Secukinumab', 'Guselkumab, Risankizumab, Tildrakizumab', 'Spesolimab',
        'Adalimumab, Etanercept, Infliximab', 'Anakinra, Canakinumab', 'Tocilizumab'])

    df['col_coord'] = df['disease'].replace({
        'Atopic dermatitis': 3, 'Bullous pemphigoid': 4, 'Granuloma annulare': 9, 'Lichen planus': 0,
        'Lupus erythematosus': 1, 'Nummular eczema': 2, 'Pityriasis rubra pilaris': 6, 'Psoriasis': 5,
        'Psoriasis pustulosa': 10, 'Pyoderma gangrenosum': 11, 'Sarcoidosis': 8, 'Scleroderma': 7})

    # Sort by pattern 1-5
    df['disease'] = df['disease'].cat.reorder_categories([
        'Lichen planus', 'Lupus erythematosus', 'Nummular eczema', 'Atopic dermatitis',
        'Bullous pemphigoid', 'Psoriasis', 'Pityriasis rubra pilaris', 'Scleroderma',
        'Sarcoidosis', 'Granuloma annulare', 'Psoriasis pustulosa', 'Pyoderma gangrenosum'])

    # Reorder categories in efficacy
    df['efficacy'] = df['efficacy'].cat.reorder_categories(
        ['very high', 'high', 'moderate', 'low', 'ineffective', 'conflict'])

    # for secondary y axis add empty string
    second_yaxis_labels = df['drug'].cat.categories  # [::-1]

    # How to show if clinical trial is still ongoing?
    minor_ticks = np.arange(0.5, len(df['disease'].unique()) - 1, 1)
    minor_yticks = np.arange(0.5, len(df['target'].unique()) - 1, 1)
    major_ticks = np.arange(0, len(df['disease'].unique()), 1)
    major_yticks = np.arange(0, len(df['target'].unique()), 1)

    df.to_excel(os.path.join(save_folder, 'Figure_1N_Table.xlsx'))

    # Plot
    fig = sns.relplot(x="disease", y="target", hue="efficacy", size='level_of_evidence', style='efficacy_markers',
                      sizes=(5, 250), alpha=1, palette=color_efficacy, height=5, data=df, markers=markers, zorder=3,
                      aspect=1.1, edgecolor="black")
    ax = fig.axes[0, 0]
    ax.set(xlim=[-0.5, 11.5], ylim=[-0.5, 10.5], autoscale_on=False)

    # Color areas of patterns
    ax.axvspan(-0.5, 1.5, facecolor='darkorange', alpha=1, zorder=1)
    ax.axvspan(1.5, 3.5, facecolor='darkred', alpha=1, zorder=1)
    ax.axvspan(3.5, 4.5, facecolor='tomato', alpha=1, zorder=1)
    ax.axvspan(4.5, 6.5, facecolor='royalblue', alpha=1, zorder=1)
    ax.axvspan(6.5, 7.5, facecolor='mediumvioletred', alpha=1, zorder=1)
    ax.axvspan(7.5, 9.5, facecolor='magenta', alpha=1, zorder=1)
    ax.axvspan(9.5, 11.5, facecolor='darkviolet', alpha=1, zorder=1)
    # Cover up with horizontal white stripes TODO
    for ind in np.arange(5, 10.5, 0.5)[::-1]:
        ax.axhspan(-0.5, ind, facecolor='white', alpha=0.2, zorder=1)

    # Add second y axis - showing the drug annotations
    par2 = ax.twinx()
    par2.spines["left"].set_position(("axes", -0.2))
    make_patch_spines_invisible(par2)
    par2.spines["left"].set_visible(True)
    par2.yaxis.set_label_position('left')
    par2.yaxis.set_ticks_position('left')
    par2.set_ylabel("drug")
    par2.set_yticks(np.linspace(0.05, 0.96, len(second_yaxis_labels)))
    par2.set_yticklabels(second_yaxis_labels)
    tkw = dict(size=4, width=1.5)
    par2.tick_params(axis='y', **tkw)

    # TODO Add colorbar instead of legend for color
    # cax = fig.fig.add_axes([ax.get_position().x1 + 0.2, ax.get_position().y0 + 0.5, 0.06, ax.get_position().height / 2])
    # ax.figure.colorbar(cmap, cax=cax)
    # Manage ticks
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True)
    ax.grid(axis='x', color='black', linestyle='-', linewidth=0.5, which='minor', alpha=0.2)
    for ind in range(len(df)):
        # print(x_pos[ind], y_pos[ind])
        if df['efficacy_combined'][ind] != '':
            ax.text(y=df['row_coord'][ind] + 0.2, x=df['col_coord'][ind] - 0.1,
                    s=df['efficacy_combined'][ind], horizontalalignment='left')
    ax.set_xticklabels(df['disease'].cat.categories, rotation=315, ha='left')
    # ax.xticks(rotation=315, ha='left')
    plt.ylabel('')
    plt.xlabel('')
    # dont show minor ticks in plot
    ax.tick_params(axis='x', which='minor', colors='white')
    ax.tick_params(axis='y', which='minor', colors='white')
    # Put a legend to the right of the current axis
    # sns.move_legend(ax, loc='center left', bbox_to_anchor=(0.6, 0.5))

    # Add second x-axis at the top of plot showing the pattern
    new_tick_locations = np.array([0.09, 0.25, 0.375, 0.5, 0.625, 0.75, 0.92])
    par_x = ax.twiny()
    par_x.spines["top"].set_position(("axes", 1.02))
    make_patch_spines_invisible(par_x)
    par_x.spines["top"].set_visible(True)
    par_x.xaxis.set_label_position('top')
    par_x.xaxis.set_ticks_position('top')
    par_x.set_xlabel("Pattern")
    par_x.set_xticks(new_tick_locations)
    par_x.set_xticklabels(['1', '2a', '2b', '3', '4a', '4b', '5'])
    tkw = dict(size=4, width=1.5)
    par_x.tick_params(axis='x', **tkw)

    # fig.tight_layout()
    # Save the figure
    fig.savefig(os.path.join(save_folder, "Figure_1N_Table_plot.pdf"), bbox_inches="tight")
    plt.close('all')


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_1N', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    input_folder = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                                'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'Paper', 'Figures',
                                'Figure_1', 'Table')

    main(input_dir=input_folder, save_folder=output_dir)


