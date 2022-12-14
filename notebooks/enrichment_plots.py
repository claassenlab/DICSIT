"""
Copyright 2022 Jan T. Schleicher
"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
from scipy.stats import mannwhitneyu


def enrichment_per_sample(meta_data: pd.DataFrame, class_col: str, sample_col: str,
                          cell_type_col: str, ratio_max_response=.5, logy=False, color_dict=None,
                          save=False, figsize=(7, 6), out_dir=""):
    """
    Plot per sample selected cell type enrichments for each filter
    @param meta_data: pandas DataFrame with cell meta data and filter response
    @param class_col: column in meta_data containing the classification classes
    @param sample_col: column in meta_data containing the sample names
    @param cell_type_col: column in meta_data containing the cell types
    @param ratio_max_response: threshold for selected cells as ratio of maximal filter response
    @param logy: log scale the y axis
    @param color_dict: dictionary with colors for each cell type
    @param save: save the plot
    @param figsize: tuple; size of the figure
    @param out_dir: path to output directory for the plots
    @return:
    """
    meta_data_df = meta_data.copy()
    filter_response_cols = meta_data.filter(regex="response_filter").columns
    
    # get selected cells
    for i, column in enumerate(filter_response_cols):
        meta_data_df[f"selected_filter_{i}"] = meta_data_df[column] >= ratio_max_response * max(meta_data_df[column])
    
    # get per sample statistics of selected cells (cell types)
    for i, column in enumerate(filter_response_cols):
        print(f"Selected cell type enrichments for filter {i}")
        if class_col is not None:
            # get class with higher number of selected cells
            per_group_selected = meta_data_df.groupby(class_col)[f"selected_filter_{i}"].sum()
            group = per_group_selected.index[per_group_selected.argmax()]
            group_data = meta_data_df[meta_data_df[class_col] == group].drop(class_col, axis=1)
            print(f"Group {group}: \n{per_group_selected[group] / per_group_selected.sum() * 100:.2f}% of "
                  f"all selected cells")
        else:
            group_data = meta_data_df
        
        # compute frequency and enrichment per sample
        per_sample_total = group_data.groupby(sample_col).size()
        per_sample_cell_type_counts = group_data.groupby([sample_col, cell_type_col]).size().unstack(fill_value=0)
        per_sample_selected = group_data.groupby(sample_col)[f"selected_filter_{i}"].sum()
        per_sample_expected = per_sample_cell_type_counts.div(per_sample_total, axis=0)\
            .multiply(per_sample_selected, axis=0)
        per_sample_selected_counts = group_data[group_data[f"selected_filter_{i}"]].groupby([sample_col, cell_type_col])\
            .size().unstack(fill_value=0)
        per_sample_enrichment_scores = per_sample_selected_counts.div(per_sample_expected).fillna(0)
        
        enrichment_scores = per_sample_enrichment_scores.reset_index().melt(id_vars=sample_col, value_name="score")
        
        # aggregate results
        aggregated_results = per_sample_enrichment_scores.agg(["median", "mad"]).unstack(level=-1)\
            .reset_index().set_axis([cell_type_col, "results", "values"], axis=1)\
            .pivot(index=cell_type_col, columns="results", values="values").reset_index()\
            .sort_values(by=["median", "mad"], ascending=False)
        cell_types_plot = aggregated_results.head(5)[cell_type_col]
        
        # make boxplots
        fig, ax = plt.subplots(figsize=figsize)
        sns.boxplot(x=cell_type_col, y="score", showfliers=False, palette=color_dict,
                    data=enrichment_scores[enrichment_scores[cell_type_col].isin(cell_types_plot)],
                    order=cell_types_plot, ax=ax)
        if enrichment_scores[cell_type_col].str.len().max() > 3:
            ax.set_xticklabels(cell_types_plot, rotation=45, ha="right")
        ax.set_xlabel("cell type")
        ax.set_ylabel("enrichment score")
        if logy:
            ax.set_yscale("log")
        ax.axhline(y=1, color="k", ls="--")
        sns.despine()
        fig.tight_layout()
        if save:
            fig.savefig(os.path.join(out_dir, f"enrichment_per_sample_{cell_type_col.replace('.','_')}_filter_{i}.svg"),
                        dpi=200)
        plt.show()
        plt.close()

        
def frequency_per_sample(meta_data: pd.DataFrame, class_col: str, sample_col: str,
                         cell_type_col: str, ratio_max_response=.5, logy=False, color_dict=None,
                         save=False, figsize=(7, 6), out_dir=""):
    """
    Plot per sample selected cell type frequencies for each filter
    @param meta_data: pandas DataFrame with cell meta data and filter response
    @param class_col: column in meta_data containing the classification classes
    @param sample_col: column in meta_data containing the sample names
    @param cell_type_col: column in meta_data containing the cell types
    @param ratio_max_response: threshold for selected cells as ratio of maximal filter response
    @param logy: log scale the y axis
    @param color_dict: dictionary with colors for each cell type
    @param save: save the plot
    @param figsize: tuple; size of the figure
    @param out_dir: path to output directory for the plots
    @return:
    """
    meta_data_df = meta_data.copy()
    filter_response_cols = meta_data.filter(regex="response_filter").columns
    
    # get selected cells
    for i, column in enumerate(filter_response_cols):
        meta_data_df[f"selected_filter_{i}"] = meta_data_df[column] >= ratio_max_response * max(meta_data_df[column])
    
    # get per sample statistics of selected cells (cell types)
    for i, column in enumerate(filter_response_cols):
        print(f"Selected cell type frequencies for filter {i}")
        if class_col is not None:
            # get class with higher number of selected cells
            per_group_selected = meta_data_df.groupby(class_col)[f"selected_filter_{i}"].sum()
            group = per_group_selected.index[per_group_selected.argmax()]
            group_data = meta_data_df[meta_data_df[class_col] == group].drop(class_col, axis=1)
            print(f"Group {group}: \n{per_group_selected[group] / per_group_selected.sum() * 100:.2f}% of "
                  f"all selected cells")
        else:
            group_data = meta_data_df
        
        # compute frequency per sample
        per_sample_selected = group_data.groupby(sample_col)[f"selected_filter_{i}"].sum()
        per_sample_selected_counts = group_data[group_data[f"selected_filter_{i}"]].groupby([sample_col, cell_type_col])\
            .size().unstack(fill_value=0)
        per_sample_frequencies = per_sample_selected_counts.div(per_sample_selected, axis=0).fillna(0)
        
        frequencies = per_sample_frequencies.reset_index().melt(id_vars=sample_col, value_name="score")
        
        # aggregate results
        aggregated_results = per_sample_frequencies.agg(["median", "mad"]).unstack(level=-1)\
            .reset_index().set_axis([cell_type_col, "results", "values"], axis=1)\
            .pivot(index=cell_type_col, columns="results", values="values").reset_index()\
            .sort_values(by=["median", "mad"], ascending=False)
        cell_types_plot = aggregated_results.head(5)[cell_type_col]
        
        # make boxplots
        fig, ax = plt.subplots(figsize=figsize)
        sns.boxplot(x=cell_type_col, y="score", showfliers=False, palette=color_dict,
                    data=frequencies[frequencies[cell_type_col].isin(cell_types_plot)],
                    order=cell_types_plot, ax=ax)
        if frequencies[cell_type_col].str.len().max() > 3:
            ax.set_xticklabels(cell_types_plot, rotation=45, ha="right")
        ax.set_xlabel("cell type")
        ax.set_ylabel("frequency")
        if logy:
            ax.set_yscale("log")
        sns.despine()
        fig.tight_layout()
        if save:
            fig.savefig(os.path.join(out_dir, f"frequency_per_sample_{cell_type_col.replace('.','_')}_filter_{i}.svg"),
                        dpi=200)
        plt.show()
        plt.close()


def plot_selected_proportion(meta_data: pd.DataFrame, cell_type_col: str, cell_type: str,
                             ratio_max_response=.5, filter_response_col="response_filter_0", 
                             palette="Set2", figsize=(3, 3), save="", out_dir=""):
    """
    Plot the proportion of a cell type that is selected
    @param meta_data: pandas DataFrame with cell meta data and filter response
    @param cell_type_col: column in meta_data containing the cell types
    @param cell_type: cell type for which proportions are plotted
    @param ratio_max_response: threshold for selected cells as ratio of maximal filter response
    @param filter_response_col: column with the filter response
    @param palette: matplotlib color palette
    @param figsize: tuple; size of the figure
    @param save: file name of output plot (w/o file ending)
    @return:
    """
    df = meta_data.copy()
    df["selected"] = df[filter_response_col] >= ratio_max_response * max(df[filter_response_col])
    df = df[df[cell_type_col] == cell_type]
    n_s = df["selected"].sum()
    n_u = len(df) - n_s
    
    colors = sns.color_palette(palette)[:3]
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.pie([n_u, n_s], explode=[0, .1], colors=colors, startangle=90)
    fig.tight_layout()
    if save:
        fig.savefig(f"{os.path.join(out_dir, save)}.svg", dpi=200)
    plt.show()
    plt.close()


def plot_selected_cell_frequencies_classification(selected_cells: pd.DataFrame, class_col: str, class_order: list,
                                                  filter_idx: int, save=False, out_dir="", palette="Set2", figsize=(4, 4)):
    """
    Plot the selected cell frequencies
    @param selected_cells: pandas DataFrame containing 'selected_filter_{filter}_freq' column
    @param class_col: column in selected_cells containing the classes
    @param class_order: list; order of classes
    @param filter_idx: int; index of filter to plot
    @param save: save the plot
    @param out_dir: path to output directory for the plots
    @param palette: matplotlib color palette
    @param figsize: tuple; size of the figure
    @return:
    """
    n_classes = len(class_order)

    if n_classes == 2:
        # compute p-value for selected cell frequencies
        _, pval = mannwhitneyu(
            selected_cells[selected_cells[class_col] == class_order[0]][f"selected_filter_{filter_idx}_freq"],
            selected_cells[selected_cells[class_col] == class_order[1]][f"selected_filter_{filter_idx}_freq"])

    # make boxplot with swarmplot
    fig, ax = plt.subplots(figsize=figsize)

    sns.boxplot(x=class_col, y=f"selected_filter_{filter_idx}_freq", data=selected_cells,
                order=class_order, palette=palette, ax=ax, width=.5, showfliers=False)
    sns.swarmplot(x=class_col, y=f"selected_filter_{filter_idx}_freq", data=selected_cells,
                  order=class_order, ax=ax, color=".25")
    ax.set_ylabel("selected population frequency [%]")

    if n_classes == 2:
        y = selected_cells[f"selected_filter_{filter_idx}_freq"].max() + 5
        h = 2

        ax.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1.2, c="k")
        ax.text(0.5, y + 2 * h, f"p = {pval:.2g}", ha="center", va="bottom")
    sns.despine()
    plt.tight_layout()
    plt.show()
    if save:
        fig.savefig(os.path.join(out_dir, f"selected_population_frequencies_filter_{filter_idx}.svg"))
