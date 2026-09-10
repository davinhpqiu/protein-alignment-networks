"""Build report macros and compact figures from authoritative saved results."""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from protein_alignment_networks import build_similarity_graph

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"
FIGURES = REPORT / "figures"


def configure_plot_style() -> None:
    """Apply one print-oriented type and colour system to report figures."""

    plt.switch_backend("Agg")
    plt.rcParams.update(
        {
            "font.size": 9.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8,
            "lines.linewidth": 1.5,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def copy_static_assets() -> None:
    sources = {
        ROOT / "Thesis_Template_OxfordPav/Logos/logo2.png": FIGURES / "logo2.png",
        ROOT / "references/references.bib": REPORT / "references.bib",
    }
    for source, destination in sources.items():
        if not source.exists():
            raise FileNotFoundError(f"required report asset is absent: {source}")
        shutil.copyfile(source, destination)


def selected_row(table: pd.DataFrame, **conditions: object) -> pd.Series:
    mask = pd.Series(True, index=table.index)
    for column, value in conditions.items():
        mask &= table[column] == value
    rows = table.loc[mask]
    if len(rows) != 1:
        raise ValueError(f"expected one row for {conditions}, found {len(rows)}")
    return rows.iloc[0]


def median_result(
    table: pd.DataFrame,
    *,
    reference: str,
    method: str,
    rule: str,
    threshold: float,
    metric: str,
) -> float:
    rows = table.loc[
        (table["reference_label"] == reference)
        & (table["method"] == method)
        & (table["threshold_rule"] == rule)
        & (table["threshold"] == threshold)
    ]
    if rows.empty:
        raise ValueError("requested report result is absent")
    return float(rows[metric].median())


def write_macros() -> None:
    summary_manifest = json.loads(
        (
            ROOT / "outputs/tables/pfam_large_panel_graph_summary/summary_manifest.json"
        ).read_text()
    )
    dataset_manifest = json.loads(
        (
            ROOT / "data/processed/pfam_large_panel/pfam_large_panel_manifest.json"
        ).read_text()
    )
    experiment_manifest = json.loads(
        (
            ROOT / "data/processed/pfam_large_panel/resampling/experiment_manifest.json"
        ).read_text()
    )
    score_manifest = json.loads(
        (ROOT / "outputs/tables/pfam_large_panel_scores/run_manifest.json").read_text()
    )
    rule = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_graph_summary/rule_summary.tsv",
        sep="\t",
    )
    paired = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_summary/paired_method_summary.tsv",
        sep="\t",
    )
    factor = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_graph_summary/factor_summary.tsv",
        sep="\t",
    )
    null = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_score_permutation_null/graph_recovery_results.tsv",
        sep="\t",
    )
    unweighted = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_unweighted_reference/graph_recovery_results.tsv",
        sep="\t",
    )
    full = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_full/graph_recovery_results.tsv",
        sep="\t",
    )

    family = "family_accession"
    bio_top5 = selected_row(
        rule,
        reference_label=family,
        method="biopython",
        threshold_rule="top_k",
        threshold=5.0,
    )
    blast_top5 = selected_row(
        rule,
        reference_label=family,
        method="blast",
        threshold_rule="top_k",
        threshold=5.0,
    )
    bio_top5_clan = selected_row(
        rule,
        reference_label="clan_accession",
        method="biopython",
        threshold_rule="top_k",
        threshold=5.0,
    )
    blast_top5_clan = selected_row(
        rule,
        reference_label="clan_accession",
        method="blast",
        threshold_rule="top_k",
        threshold=5.0,
    )
    bio_p90 = selected_row(
        rule,
        reference_label=family,
        method="biopython",
        threshold_rule="percentile",
        threshold=0.9,
    )
    blast_p90 = selected_row(
        rule,
        reference_label=family,
        method="blast",
        threshold_rule="percentile",
        threshold=0.9,
    )
    matched_two = selected_row(
        paired,
        reference_label=family,
        threshold_rule="target_density",
        threshold=0.02,
    )
    observed_top5_ami = median_result(
        full.loc[
            (full["clans_per_collection"] == 4)
            & (full["families_per_clan"] == 2)
            & (full["domains_per_family"] == 20)
        ],
        reference=family,
        method="biopython",
        rule="top_k",
        threshold=5.0,
        metric="ami",
    )
    reference_full = full.loc[
        (full["clans_per_collection"] == 4)
        & (full["families_per_clan"] == 2)
        & (full["domains_per_family"] == 20)
    ]
    blast_p99_nmi = median_result(
        reference_full,
        reference=family,
        method="blast",
        rule="percentile",
        threshold=0.99,
        metric="nmi",
    )
    blast_p99_ami = median_result(
        reference_full,
        reference=family,
        method="blast",
        rule="percentile",
        threshold=0.99,
        metric="ami",
    )
    blast_p99_communities = median_result(
        reference_full,
        reference=family,
        method="blast",
        rule="percentile",
        threshold=0.99,
        metric="inferred_community_count",
    )
    null_top5_ami = median_result(
        null,
        reference=family,
        method="biopython",
        rule="top_k",
        threshold=5.0,
        metric="ami",
    )
    unweighted_top5_ami = median_result(
        unweighted,
        reference=family,
        method="biopython",
        rule="top_k",
        threshold=5.0,
        metric="ami",
    )
    blast_audit = score_manifest["parameters"]["blast"]["directionality_audit"]
    composition = factor.loc[
        (factor["reference_label"] == family)
        & (factor["method"] == "biopython")
        & (factor["threshold_rule"] == "top_k")
        & (factor["threshold"] == 5)
    ]
    blast_composition = factor.loc[
        (factor["reference_label"] == family)
        & (factor["method"] == "blast")
        & (factor["threshold_rule"] == "top_k")
        & (factor["threshold"] == 5)
    ]
    representative = full.loc[
        (full["collection_id"] == "c4_f2_d20_r001")
        & (full["reference_label"] == family)
        & (full["threshold_rule"] == "top_k")
        & (full["threshold"] == 5)
    ].set_index("method")
    condition_count = (
        len(experiment_manifest["clans_per_collection"])
        * len(experiment_manifest["families_per_clan"])
        * len(experiment_manifest["domains_per_family"])
    )
    replicates = experiment_manifest["replicates_per_condition"]
    macro_values = {
        "MasterDomains": dataset_manifest["selected_sequence_count"],
        "MasterFamilies": len(dataset_manifest["families"]),
        "MasterClans": len(
            {item["clan_accession"] for item in dataset_manifest["families"]}
        ),
        "MasterPairs": dataset_manifest["unique_pair_count"],
        "SampledCollections": experiment_manifest["collection_count"],
        "GraphsEvaluated": summary_manifest["graph_count"],
        "PrimaryCollections": summary_manifest["primary_analysis_collection_count"],
        "ExcludedCells": summary_manifest["excluded_composition_cell_count"],
        "ExcludedCollections": (
            experiment_manifest["collection_count"]
            - summary_manifest["primary_analysis_collection_count"]
        ),
        "PrimaryConditions": (
            condition_count - summary_manifest["excluded_composition_cell_count"]
        ),
        "ConditionCount": condition_count,
        "ReplicatesPerCondition": replicates,
        "SamplingSeed": experiment_manifest["base_seed"],
        "TargetPerFamily": dataset_manifest["target_per_family"],
        "BioTopFiveAMI": f"{bio_top5['median_ami']:.3f}",
        "BlastTopFiveAMI": f"{blast_top5['median_ami']:.3f}",
        "BioTopFiveClanAMI": f"{bio_top5_clan['median_ami']:.3f}",
        "BlastTopFiveClanAMI": f"{blast_top5_clan['median_ami']:.3f}",
        "BioTopFiveDensityPct": f"{100 * bio_top5['median_density']:.2f}",
        "BlastTopFiveDensityPct": f"{100 * blast_top5['median_density']:.2f}",
        "BioPercentileNinetyDensityPct": f"{100 * bio_p90['median_density']:.2f}",
        "BlastPercentileNinetyDensityPct": f"{100 * blast_p90['median_density']:.2f}",
        "MatchedTwoAbsDeltaNMI": f"{matched_two['median_absolute_nmi_difference']:.4f}",
        "ObservedTopFiveAMI": f"{observed_top5_ami:.3f}",
        "NullTopFiveAMI": f"{null_top5_ami:.3f}",
        "UnweightedTopFiveAMI": f"{unweighted_top5_ami:.3f}",
        "BlastPercentileNinetyNineNMI": f"{blast_p99_nmi:.3f}",
        "BlastPercentileNinetyNineAMI": f"{blast_p99_ami:.3f}",
        "BlastPercentileNinetyNineCommunities": (f"{blast_p99_communities:.0f}"),
        "BlastReportedPairs": blast_audit["undirected_pair_count"],
        "BlastCoveragePct": (
            f"{100 * blast_audit['undirected_pair_count'] / dataset_manifest['unique_pair_count']:.2f}"
        ),
        "CompositionMinAMI": f"{composition['median_ami'].min():.3f}",
        "CompositionMaxAMI": f"{composition['median_ami'].max():.3f}",
        "BlastCompositionMinAMI": f"{blast_composition['median_ami'].min():.3f}",
        "BlastCompositionMaxAMI": f"{blast_composition['median_ami'].max():.3f}",
        "RepresentativeNodes": int(representative.loc["biopython", "nodes"]),
        "RepresentativeFamilies": int(
            representative.loc["biopython", "reference_group_count"]
        ),
        "RepresentativeBioEdges": int(representative.loc["biopython", "edges"]),
        "RepresentativeBlastEdges": int(representative.loc["blast", "edges"]),
        "RepresentativeBioAMI": (f"{representative.loc['biopython', 'ami']:.3f}"),
        "RepresentativeBlastAMI": f"{representative.loc['blast', 'ami']:.3f}",
    }
    lines = [
        "% Generated by scripts/build_report_assets.py; do not edit by hand.",
        "% Sources: dataset, experiment, score, graph-summary, null, and unweighted manifests/TSVs.",
    ]
    lines.extend(
        f"\\newcommand{{\\{name}}}{{{value}}}" for name, value in macro_values.items()
    )
    (REPORT / "generated_results.tex").write_text("\n".join(lines) + "\n")


def plot_percentile_mismatch() -> None:
    rule = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_graph_summary/rule_summary.tsv",
        sep="\t",
    )
    percentile = rule.loc[rule["threshold_rule"] == "percentile"].copy()
    method_labels = {"biopython": "Biopython", "blast": "BLAST"}
    colours = {"biopython": "#002147", "blast": "#B24A35"}
    figure, axes = plt.subplots(1, 3, figsize=(10.2, 3.2))
    panels = [
        ("family_accession", "median_density", "Realised density", "Density"),
        ("family_accession", "median_nmi", "Family recovery", "Median NMI"),
        ("clan_accession", "median_nmi", "Clan recovery", "Median NMI"),
    ]
    for axis, (reference, metric, title, ylabel) in zip(axes, panels, strict=True):
        panel = percentile.loc[percentile["reference_label"] == reference]
        for method, group in panel.groupby("method"):
            ordered = group.sort_values("threshold", ascending=False)
            values = ordered[metric]
            if metric == "median_density":
                values = 100 * values
            axis.plot(
                100 * (1 - ordered["threshold"]),
                values,
                color=colours[method],
                marker="o",
                linewidth=1.5,
                label=method_labels[method],
            )
        axis.set_title(title)
        axis.set_xlabel("Available scores retained (%)")
        axis.set_ylabel(ylabel + (" (%)" if metric == "median_density" else ""))
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(
        FIGURES / "percentile_density_mismatch.png",
        dpi=240,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_rule_recovery() -> None:
    rule = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_graph_summary/rule_summary.tsv",
        sep="\t",
    )
    subset = rule.loc[rule["reference_label"] == "family_accession"].copy()
    colours = {
        "top_k": "#002147",
        "target_density": "#B24A35",
        "percentile": "#6B8E23",
    }
    markers = {"biopython": "o", "blast": "s"}
    method_labels = {"biopython": "Biopython", "blast": "BLAST"}
    rule_labels = {
        "top_k": r"top-$k$",
        "target_density": "target density",
        "percentile": "percentile",
    }
    figure, axis = plt.subplots(figsize=(7.2, 4.4))
    for (method, graph_rule), group in subset.groupby(["method", "threshold_rule"]):
        ordered = group.sort_values("median_density")
        axis.plot(
            100 * ordered["median_density"],
            ordered["median_ami"],
            color=colours[graph_rule],
            marker=markers[method],
            linewidth=1.4,
            markersize=5,
            label=f"{method_labels[method]} — {rule_labels[graph_rule]}",
        )
    axis.set_xlabel("Median realised edge density (%)")
    axis.set_ylabel("Median family AMI")
    axis.set_ylim(0, 1)
    axis.grid(alpha=0.2)
    axis.legend(ncol=2, fontsize=8, frameon=False)
    figure.tight_layout()
    figure.savefig(FIGURES / "rule_recovery.png", dpi=240, bbox_inches="tight")
    plt.close(figure)


def _packed_component_layout(
    graphs: tuple[nx.Graph, ...],
    *,
    seed: int,
) -> dict[str, np.ndarray]:
    """Pack union-graph components, then lay out nodes within each component."""

    union = nx.Graph()
    for graph in graphs:
        union.add_nodes_from(graph.nodes)
        union.add_edges_from(graph.edges)
    components = sorted(
        nx.connected_components(union),
        key=lambda nodes: (-len(nodes), min(map(str, nodes))),
    )
    columns = max(1, math.ceil(math.sqrt(2 * len(components))))
    rows = math.ceil(len(components) / columns)
    positions: dict[str, np.ndarray] = {}
    for index, nodes in enumerate(components):
        subgraph = union.subgraph(nodes)
        if len(nodes) == 1:
            local = {next(iter(nodes)): np.zeros(2)}
        else:
            local = nx.spring_layout(
                subgraph,
                seed=seed + index,
                iterations=300,
                weight=None,
                scale=0.78,
            )
        column = index % columns
        row = index // columns
        centre = np.array(
            [2.0 * (column - (columns - 1) / 2), 2.0 * ((rows - 1) / 2 - row)]
        )
        positions.update(
            {node: coordinate + centre for node, coordinate in local.items()}
        )
    return positions


def plot_representative_graphs() -> None:
    collection_id = "c4_f2_d20_r001"
    membership = pd.read_csv(
        ROOT / "data/processed/pfam_large_panel/resampling/collection_membership.tsv",
        sep="\t",
    )
    sample = membership.loc[membership["collection_id"] == collection_id].copy()
    if sample.empty:
        raise ValueError(f"representative collection is absent: {collection_id}")
    protein_ids = sample["node_id"].tolist()
    pairs = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_scores/paired_scores.tsv",
        sep="\t",
    )
    in_collection = pairs["protein_a"].isin(protein_ids) & pairs["protein_b"].isin(
        protein_ids
    )
    pairs = pairs.loc[in_collection].copy()
    metadata = sample[["node_id", "family_accession", "clan_accession"]]
    graphs = {
        method: build_similarity_graph(
            pairs,
            protein_ids,
            method=method,
            threshold=5,
            threshold_rule="top_k",
            node_metadata=metadata,
            node_id_column="node_id",
        )
        for method in ("biopython", "blast")
    }
    recovery = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_full/graph_recovery_results.tsv",
        sep="\t",
    )
    expected = recovery.loc[
        (recovery["collection_id"] == collection_id)
        & (recovery["reference_label"] == "family_accession")
        & (recovery["threshold_rule"] == "top_k")
        & (recovery["threshold"] == 5)
    ].set_index("method")
    for method, graph in graphs.items():
        if graph.number_of_edges() != int(expected.loc[method, "edges"]):
            raise ValueError(
                f"{method} representative graph does not match saved result"
            )

    positions = _packed_component_layout(tuple(graphs.values()), seed=42)
    families = sorted(sample["family_accession"].unique())
    palette = plt.get_cmap("tab10")
    colours = {family: palette(index) for index, family in enumerate(families)}
    figure, axes = plt.subplots(2, 1, figsize=(8.4, 8.2), sharex=True, sharey=True)
    labels = {"biopython": "Biopython", "blast": "BLAST"}
    for axis, (method, graph) in zip(axes, graphs.items(), strict=True):
        nx.draw_networkx_edges(
            graph,
            positions,
            ax=axis,
            edge_color="#65737E",
            width=0.75,
            alpha=0.58,
        )
        nx.draw_networkx_nodes(
            graph,
            positions,
            ax=axis,
            node_color=[
                colours[str(graph.nodes[node]["family_accession"])] for node in graph
            ],
            node_size=24,
            linewidths=0.3,
            edgecolors="white",
        )
        axis.set_title(
            f"{labels[method]}: {graph.number_of_edges()} edges; "
            f"family AMI {expected.loc[method, 'ami']:.3f}",
        )
        axis.set_aspect("equal", adjustable="box")
        axis.set_facecolor("#FCFCFC")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_color("#D6DBE1")
            spine.set_linewidth(0.7)
    legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=colours[family],
            markeredgecolor="none",
            markersize=6,
            label=family,
        )
        for family in families
    ]
    figure.legend(
        handles=legend,
        loc="lower center",
        ncol=4,
        frameon=False,
        fontsize=8,
        title="Pfam family",
        title_fontsize=8.5,
    )
    figure.tight_layout(rect=(0, 0.12, 1, 1), h_pad=0.8)
    figure.savefig(
        FIGURES / "representative_graphs.png",
        dpi=260,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_matched_density_differences() -> None:
    paired = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_summary/paired_method_differences.tsv",
        sep="\t",
    )
    paired = paired.loc[
        paired["comparison_status"] == "primary_matched_edge_budget"
    ].copy()
    thresholds = [0.005, 0.01, 0.02]
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
    references = [
        ("family_accession", "Family labels"),
        ("clan_accession", "Clan labels"),
    ]
    for axis, (reference, title) in zip(axes, references, strict=True):
        panel = paired.loc[paired["reference_label"] == reference]
        values = [
            panel.loc[panel["threshold"] == threshold, "nmi_blast_minus_biopython"]
            for threshold in thresholds
        ]
        boxplot = axis.boxplot(
            values,
            tick_labels=["0.5%", "1%", "2%"],
            patch_artist=True,
            showfliers=False,
            widths=0.48,
        )
        for box in boxplot["boxes"]:
            box.set_facecolor("#DCE6F1")
            box.set_edgecolor("#002147")
        for median in boxplot["medians"]:
            median.set_color("#B24A35")
            median.set_linewidth(1.5)
        axis.axhline(0, color="#333333", linewidth=0.8)
        axis.set_title(title)
        axis.set_xlabel("Target edge density")
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("BLAST NMI minus Biopython NMI")
    figure.tight_layout()
    figure.savefig(
        FIGURES / "matched_density_differences.png",
        dpi=240,
        bbox_inches="tight",
    )
    plt.close(figure)


def _rule_order() -> list[tuple[str, float, str]]:
    return [
        ("percentile", 0.90, "p90"),
        ("percentile", 0.95, "p95"),
        ("percentile", 0.98, "p98"),
        ("percentile", 0.99, "p99"),
        ("target_density", 0.005, "d0.5%"),
        ("target_density", 0.01, "d1%"),
        ("target_density", 0.02, "d2%"),
        ("target_density", 0.05, "d5%"),
        ("top_k", 1.0, "k1"),
        ("top_k", 2.0, "k2"),
        ("top_k", 5.0, "k5"),
        ("top_k", 10.0, "k10"),
    ]


def plot_rule_frequency(metric: str) -> None:
    """Plot one agreement distribution for all rules at reference composition."""

    if metric not in {"nmi", "ami"}:
        raise ValueError("metric must be 'nmi' or 'ami'")
    full = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_full/graph_recovery_results.tsv",
        sep="\t",
    )
    matched = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_matched_density/graph_recovery_results.tsv",
        sep="\t",
    )
    results = pd.concat([full, matched], ignore_index=True)
    results = results.loc[
        (results["clans_per_collection"] == 4)
        & (results["families_per_clan"] == 2)
        & (results["domains_per_family"] == 20)
    ]
    order = _rule_order()
    colours = plt.colormaps["tab20"](np.linspace(0, 1, len(order)))
    linestyles = {"percentile": "-", "target_density": "--", "top_k": ":"}
    references = [
        ("family_accession", "Family labels"),
        ("clan_accession", "Clan labels"),
    ]
    methods = [("biopython", "Biopython"), ("blast", "BLAST")]
    bins = np.linspace(0, 1, 26)
    centres = (bins[:-1] + bins[1:]) / 2
    figure, axes = plt.subplots(2, 2, figsize=(8.4, 6.8), sharex=True, sharey=True)
    for row, (reference, reference_title) in enumerate(references):
        for column, (method, method_title) in enumerate(methods):
            axis = axes[row, column]
            panel = results.loc[
                (results["reference_label"] == reference)
                & (results["method"] == method)
            ]
            for colour, (rule, threshold, label) in zip(colours, order, strict=True):
                values = panel.loc[
                    (panel["threshold_rule"] == rule)
                    & (panel["threshold"] == threshold),
                    metric,
                ]
                if len(values) != 50:
                    raise ValueError(
                        f"expected 50 {method} {reference} {label} values; "
                        f"found {len(values)}"
                    )
                counts, _ = np.histogram(values, bins=bins)
                axis.plot(
                    centres,
                    100 * counts / len(values),
                    color=colour,
                    linestyle=linestyles[rule],
                    linewidth=1.45,
                    label=label,
                )
            axis.set_title(f"{method_title}: {reference_title}")
            axis.set_xlim(0, 1)
            axis.grid(alpha=0.2)
    axes[0, 0].set_ylim(bottom=0)
    for axis in axes[-1]:
        axis.set_xlabel(metric.upper())
    for axis in axes[:, 0]:
        axis.set_ylabel("Collections in bin (%)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        frameon=False,
        fontsize=7.5,
    )
    figure.tight_layout(rect=(0, 0.12, 1, 1), h_pad=1.1)
    figure.savefig(
        FIGURES / f"rule_{metric}_frequency.png",
        dpi=260,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_rule_distributions() -> None:
    full = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_full/graph_recovery_results.tsv",
        sep="\t",
    )
    matched = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_matched_density/graph_recovery_results.tsv",
        sep="\t",
    )
    results = pd.concat([full, matched], ignore_index=True)
    results = results.loc[
        (results["clans_per_collection"] == 4)
        & (results["families_per_clan"] == 2)
        & (results["domains_per_family"] == 20)
        & (results["reference_label"] == "family_accession")
    ]
    order = _rule_order()
    colours = {
        "percentile": "#6B8E23",
        "target_density": "#B24A35",
        "top_k": "#002147",
    }
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), sharey=True)
    for axis, method in zip(axes, ["biopython", "blast"], strict=True):
        panel = results.loc[results["method"] == method]
        values = [
            panel.loc[
                (panel["threshold_rule"] == rule) & (panel["threshold"] == threshold),
                "ami",
            ]
            for rule, threshold, _ in order
        ]
        boxplot = axis.boxplot(
            values,
            tick_labels=[label for _, _, label in order],
            patch_artist=True,
            showfliers=False,
            widths=0.62,
        )
        for box, (rule, _, _) in zip(boxplot["boxes"], order, strict=True):
            box.set_facecolor(colours[rule])
            box.set_alpha(0.58)
        for median in boxplot["medians"]:
            median.set_color("#111111")
            median.set_linewidth(1.2)
        axis.set_title("Biopython" if method == "biopython" else "BLAST")
        axis.tick_params(axis="x", labelrotation=45)
        axis.grid(axis="y", alpha=0.2)
        axis.set_ylim(-0.05, 1.02)
    axes[0].set_ylabel("Family AMI")
    legend = [
        Line2D([0], [0], color=colour, linewidth=7, alpha=0.58, label=label)
        for label, colour in [
            ("percentile", colours["percentile"]),
            ("target density", colours["target_density"]),
            (r"top-$k$", colours["top_k"]),
        ]
    ]
    figure.legend(handles=legend, loc="lower center", ncol=3, frameon=False)
    figure.tight_layout(rect=(0, 0.08, 1, 1))
    figure.savefig(
        FIGURES / "rule_distributions.png",
        dpi=240,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_metric_comparison() -> None:
    factor = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_graph_summary/factor_summary.tsv",
        sep="\t",
    )
    subset = factor.loc[
        (factor["clans_per_collection"] == 4)
        & (factor["families_per_clan"] == 2)
        & (factor["domains_per_family"] == 20)
    ].copy()
    order = _rule_order()
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), sharex=True, sharey=True)
    references = [
        ("family_accession", "Family labels"),
        ("clan_accession", "Clan labels"),
    ]
    methods = [("biopython", "Biopython"), ("blast", "BLAST")]
    for row, (reference, reference_label) in enumerate(references):
        for column, (method, method_label) in enumerate(methods):
            axis = axes[row, column]
            panel = subset.loc[
                (subset["reference_label"] == reference) & (subset["method"] == method)
            ]
            ordered = pd.DataFrame(
                [
                    selected_row(
                        panel,
                        threshold_rule=rule,
                        threshold=threshold,
                    )
                    for rule, threshold, _ in order
                ]
            )
            x = np.arange(len(order))
            axis.plot(x, ordered["median_nmi"], marker="o", label="NMI")
            axis.plot(x, ordered["median_ami"], marker="o", label="AMI")
            axis.fill_between(
                x,
                ordered["median_nmi"].to_numpy(dtype=float),
                ordered["median_ami"].to_numpy(dtype=float),
                color="#DCE6F1",
                alpha=0.65,
            )
            axis.set_title(f"{method_label}: {reference_label}")
            axis.grid(alpha=0.2)
            axis.set_ylim(-0.05, 1.02)
            axis.set_xticks(x, [label for _, _, label in order], rotation=45)
    axes[0, 0].set_ylabel("Median agreement")
    axes[1, 0].set_ylabel("Median agreement")
    axes[0, 0].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(
        FIGURES / "metric_comparison.png",
        dpi=240,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_sensitivity_checks() -> None:
    null = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_score_permutation_null/graph_recovery_results.tsv",
        sep="\t",
    )
    unweighted = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_unweighted_reference/graph_recovery_results.tsv",
        sep="\t",
    )
    full = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_full/graph_recovery_results.tsv",
        sep="\t",
    )
    matched = pd.read_csv(
        ROOT
        / "outputs/tables/pfam_large_panel_graph_recovery_matched_density/graph_recovery_results.tsv",
        sep="\t",
    )
    family = "family_accession"
    all_recovery = pd.concat([full, matched], ignore_index=True)
    real = all_recovery
    real = real.loc[
        (real["collection_id"].isin(null["collection_id"]))
        & (real["reference_label"] == family)
        & (real["method"] == "biopython")
        & (
            ((real["threshold_rule"] == "top_k") & (real["threshold"] == 5))
            | (
                (real["threshold_rule"] == "target_density")
                & (real["threshold"] == 0.02)
            )
        )
    ].copy()
    null_family = null.loc[null["reference_label"] == family].copy()
    real_summary = real.groupby(["threshold_rule", "threshold"])["ami"].median()
    null_summary = null_family.groupby(["threshold_rule", "threshold"])["ami"].median()

    reference_real = all_recovery.loc[
        (all_recovery["clans_per_collection"] == 4)
        & (all_recovery["families_per_clan"] == 2)
        & (all_recovery["domains_per_family"] == 20)
        & (all_recovery["reference_label"] == family)
        & (
            (
                (all_recovery["threshold_rule"] == "top_k")
                & (all_recovery["threshold"] == 5)
            )
            | (
                (all_recovery["threshold_rule"] == "target_density")
                & (all_recovery["threshold"] == 0.02)
            )
        )
    ]
    weighted_summary = reference_real.groupby(
        ["method", "threshold_rule", "threshold"]
    )["ami"].median()
    unweighted_summary = (
        unweighted.loc[unweighted["reference_label"] == family]
        .groupby(["method", "threshold_rule", "threshold"])["ami"]
        .median()
    )

    figure, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    settings = [("top_k", 5.0), ("target_density", 0.02)]
    labels = ["top-k 5", "density 2%"]
    x = np.arange(len(settings))
    axes[0].bar(
        x - 0.18,
        [real_summary.loc[item] for item in settings],
        width=0.36,
        label="observed scores",
        color="#002147",
    )
    axes[0].bar(
        x + 0.18,
        [null_summary.loc[item] for item in settings],
        width=0.36,
        label="permuted scores",
        color="#B24A35",
    )
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Median family AMI")
    axes[0].set_title("Score-permutation control")
    axes[0].legend(frameon=False, fontsize=8)

    categories = [
        ("biopython", "top_k", 5.0),
        ("blast", "top_k", 5.0),
        ("biopython", "target_density", 0.02),
        ("blast", "target_density", 0.02),
    ]
    category_labels = [
        "Biopython\ntop-k",
        "BLAST\ntop-k",
        "Biopython\ndensity",
        "BLAST\ndensity",
    ]
    x = np.arange(len(categories))
    axes[1].plot(
        x,
        [weighted_summary.loc[item] for item in categories],
        marker="o",
        label="native weights",
        color="#002147",
    )
    axes[1].plot(
        x,
        [unweighted_summary.loc[item] for item in categories],
        marker="s",
        label="unweighted",
        color="#B24A35",
    )
    axes[1].set_xticks(x, category_labels)
    axes[1].set_title("Community-weight sensitivity")
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_ylim(-0.05, 1)
        axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(FIGURES / "sensitivity_checks.png", dpi=240, bbox_inches="tight")
    plt.close(figure)


def plot_composition_recovery() -> None:
    factor = pd.read_csv(
        ROOT / "outputs/tables/pfam_large_panel_graph_summary/factor_summary.tsv",
        sep="\t",
    )
    subset = factor.loc[
        (factor["reference_label"] == "family_accession")
        & (factor["threshold_rule"] == "top_k")
        & (factor["threshold"] == 5)
    ].copy()
    colours = {2: "#002147", 4: "#B24A35", 8: "#6B8E23"}
    figure, axes = plt.subplots(2, 2, figsize=(9.2, 6.2), sharex=True, sharey=True)
    methods = [("biopython", "Biopython"), ("blast", "BLAST")]
    for row, (method, method_label) in enumerate(methods):
        for column, families in enumerate([2, 3]):
            axis = axes[row, column]
            panel = subset.loc[
                (subset["method"] == method) & (subset["families_per_clan"] == families)
            ]
            for clans, group in panel.groupby("clans_per_collection"):
                ordered = group.sort_values("domains_per_family")
                axis.plot(
                    ordered["domains_per_family"],
                    ordered["median_ami"],
                    color=colours[int(clans)],
                    marker="o",
                    linewidth=1.5,
                    label=f"{int(clans)} clans",
                )
            axis.set_title(f"{method_label}: {families} families per clan")
            axis.set_xticks([10, 20, 40])
            axis.set_ylim(0.85, 1.01)
            axis.grid(alpha=0.2)
    axes[1, 0].set_xlabel("Domains per family")
    axes[1, 1].set_xlabel("Domains per family")
    axes[0, 0].set_ylabel("Median family AMI")
    axes[1, 0].set_ylabel("Median family AMI")
    axes[0, 1].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURES / "composition_recovery.png", dpi=240, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    configure_plot_style()
    copy_static_assets()
    write_macros()
    plot_representative_graphs()
    plot_percentile_mismatch()
    plot_rule_frequency("nmi")
    plot_rule_frequency("ami")
    plot_rule_recovery()
    plot_matched_density_differences()
    plot_sensitivity_checks()
    plot_composition_recovery()
    print(f"Report assets written under {REPORT}")


if __name__ == "__main__":
    main()
