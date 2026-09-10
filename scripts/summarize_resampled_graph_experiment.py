"""Create the report-facing summaries for the repeated Pfam graph experiment.

The script combines the raw top-k/percentile and matched-density result tables,
audits repeated node sets, restricts primary distributional summaries to fully
independent composition cells, calculates paired method differences, and
writes all summary TSVs and final PNG figures.  Notebook 03 loads these saved
artifacts and explains their interpretation; numerical findings should not be
manually copied into notebook Markdown.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COMPOSITION_COLUMNS = [
    "clans_per_collection",
    "families_per_clan",
    "domains_per_family",
]
REFERENCE_LABELS = (
    ("family_accession", "Family"),
    ("clan_accession", "Clan"),
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collection_fingerprints(membership: pd.DataFrame) -> pd.DataFrame:
    """Return a stable node-set fingerprint for every sampled collection."""

    required = {"collection_id", "node_id"}
    missing = required - set(membership.columns)
    if missing:
        raise ValueError(f"membership is missing: {', '.join(sorted(missing))}")

    def fingerprint(nodes: pd.Series) -> str:
        payload = "\0".join(sorted(nodes.astype(str))).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    return (
        membership.groupby("collection_id", as_index=False)
        .agg(
            node_count=("node_id", "nunique"),
            node_set_sha256=("node_id", fingerprint),
        )
        .sort_values("collection_id", ignore_index=True)
    )


def paired_method_differences(results: pd.DataFrame) -> pd.DataFrame:
    """Pair Biopython and BLAST results on the exact graph design and sample."""

    index = [
        "collection_id",
        "clans_per_collection",
        "families_per_clan",
        "domains_per_family",
        "node_count",
        "threshold_rule",
        "threshold",
        "reference_label",
    ]
    metrics = ["nmi", "ami", "density", "edges", "connected_components"]
    wide = results.pivot(index=index, columns="method", values=metrics)
    required = {
        (metric, method) for metric in metrics for method in ("biopython", "blast")
    }
    missing = required - set(wide.columns)
    if missing:
        raise ValueError("both methods are required for every paired metric")
    wide = wide.reset_index()
    wide.columns = [
        first if not second else f"{first}_{second}" for first, second in wide.columns
    ]
    for metric in metrics:
        wide[f"{metric}_blast_minus_biopython"] = (
            wide[f"{metric}_blast"] - wide[f"{metric}_biopython"]
        )
    return wide


def quantile_summary(results: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    return (
        results.groupby(groups, as_index=False)
        .agg(
            graphs=("graph_id", "nunique"),
            median_nmi=("nmi", "median"),
            lower_quartile_nmi=("nmi", lambda values: values.quantile(0.25)),
            upper_quartile_nmi=("nmi", lambda values: values.quantile(0.75)),
            median_ami=("ami", "median"),
            median_density=("density", "median"),
            median_edges=("edges", "median"),
            median_components=("connected_components", "median"),
            median_communities=("inferred_community_count", "median"),
        )
        .sort_values(groups, ignore_index=True)
    )


def rule_label(rule: str, threshold: float) -> str:
    """Return a compact, unambiguous label for a graph setting."""

    if rule == "top_k":
        return f"top-k {int(threshold)}"
    if rule == "percentile":
        return f"percentile {threshold:.2f}"
    if rule == "target_density":
        return f"density {100 * threshold:g}%"
    return f"{rule} {threshold:g}"


def method_label(method: str) -> str:
    return "BLAST" if method == "blast" else method.title()


def fully_matched_target_densities(results: pd.DataFrame) -> list[float]:
    """Return target densities supplied without an edge shortfall anywhere."""

    required = {"threshold_rule", "threshold", "selection_shortfall"}
    missing = required - set(results.columns)
    if missing:
        raise ValueError(f"results are missing: {', '.join(sorted(missing))}")
    target = results.loc[results["threshold_rule"] == "target_density"]
    if target.empty:
        return []
    shortfall = target.groupby("threshold")["selection_shortfall"].max()
    return sorted(float(value) for value in shortfall.index[shortfall == 0])


def comparison_status(
    rule: str,
    threshold: float,
    matched_target_densities: list[float],
) -> str:
    """Name the evidential role of each paired graph comparison."""

    if rule == "target_density":
        if any(np.isclose(threshold, value) for value in matched_target_densities):
            return "primary_matched_edge_budget"
        return "stress_test_edge_budget_shortfall"
    if rule == "percentile":
        return "diagnostic_unmatched_available_fraction"
    if rule == "top_k":
        return "local_neighbour_rule_realised_density"
    raise ValueError(f"unknown graph rule: {rule}")


def add_analysis_eligibility(
    results: pd.DataFrame,
    collections: pd.DataFrame,
    uniqueness: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mark composition cells suitable for repeated-sample inference.

    A cell is primary-analysis eligible only when all declared replicates have
    distinct node sets. Repeated node sets remain in the raw experiment and in
    audit figures, but do not contribute to distributional summaries.
    """

    status = uniqueness[COMPOSITION_COLUMNS + ["analysis_eligible"]]
    collection_status = collections.merge(
        status,
        on=COMPOSITION_COLUMNS,
        how="left",
        validate="many_to_one",
    )
    if collection_status["analysis_eligible"].isna().any():
        raise ValueError("analysis eligibility is missing for one or more collections")
    labelled_results = results.merge(
        collection_status[["collection_id", "analysis_eligible"]],
        on="collection_id",
        how="left",
        validate="many_to_one",
    )
    if labelled_results["analysis_eligible"].isna().any():
        raise ValueError("analysis eligibility is missing for one or more results")
    return labelled_results, collection_status


def plot_top_k_input_sensitivity(results: pd.DataFrame, output: Path) -> None:
    subset = results.loc[
        (results["method"] == "biopython") & (results["threshold_rule"] == "top_k")
    ]
    for family_count in sorted(subset["families_per_clan"].unique()):
        figure, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey="row")
        for row, (reference, label) in enumerate(
            (("family_accession", "Family NMI"), ("clan_accession", "Clan NMI"))
        ):
            for column, clan_count in enumerate((2, 4, 8)):
                axis = axes[row, column]
                panel = subset.loc[
                    (subset["reference_label"] == reference)
                    & (subset["families_per_clan"] == family_count)
                    & (subset["clans_per_collection"] == clan_count)
                ]
                medians = (
                    panel.groupby(["domains_per_family", "threshold"])["nmi"]
                    .median()
                    .unstack("threshold")
                )
                for threshold in medians.columns:
                    axis.plot(
                        medians.index,
                        medians[threshold],
                        marker="o",
                        label=f"top_k={int(threshold)}",
                    )
                axis.set_title(f"{clan_count} clans")
                axis.set_xlabel("Domains sampled per family")
                axis.set_ylabel(label)
                axis.set_ylim(0, 1.02)
                axis.set_xticks(sorted(panel["domains_per_family"].unique()))
        axes[0, 0].legend(loc="lower left")
        figure.suptitle(
            f"Biopython top-k sensitivity with {family_count} families per clan"
        )
        figure.tight_layout()
        figure.savefig(
            output / f"top_k_input_sensitivity_families_{family_count}.png",
            dpi=200,
        )
        plt.close(figure)


def plot_matched_density_method_deltas(
    paired: pd.DataFrame,
    output: Path,
    matched_target_densities: list[float],
) -> None:
    subset = paired.loc[
        (paired["threshold_rule"] == "target_density")
        & (paired["threshold"].isin(matched_target_densities))
    ]
    thresholds = sorted(subset["threshold"].unique())
    figure, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for axis, (reference, title) in zip(
        axes,
        (
            ("family_accession", "Family recovery"),
            ("clan_accession", "Clan recovery"),
        ),
    ):
        panel = subset.loc[subset["reference_label"] == reference]
        values = [
            panel.loc[
                panel["threshold"] == threshold,
                "nmi_blast_minus_biopython",
            ]
            for threshold in thresholds
        ]
        axis.boxplot(
            values, tick_labels=[f"{100 * x:g}%" for x in thresholds], showfliers=False
        )
        axis.axhline(0, color="black", linewidth=1)
        axis.set_title(title)
        axis.set_xlabel("Matched target edge density")
        axis.set_ylabel("BLAST NMI − Biopython NMI")
    figure.suptitle("Paired method differences at fully supplied target densities")
    figure.tight_layout()
    figure.savefig(output / "matched_density_method_differences.png", dpi=200)
    plt.close(figure)


def plot_percentile_mismatch(results: pd.DataFrame, output: Path) -> None:
    subset = results.loc[results["threshold_rule"] == "percentile"].copy()
    subset["retained_available_fraction"] = 1 - subset["threshold"]
    figure, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for method, group in subset.groupby("method"):
        density = group.groupby("retained_available_fraction")["density"].median()
        axes[0].plot(100 * density.index, density.values, marker="o", label=method)
        for axis, reference in zip(axes[1:], ("family_accession", "clan_accession")):
            recovery = (
                group.loc[group["reference_label"] == reference]
                .groupby("retained_available_fraction")["nmi"]
                .median()
            )
            axis.plot(100 * recovery.index, recovery.values, marker="o", label=method)
    axes[0].set_title("Realised graph density")
    axes[0].set_ylabel("Median realised density")
    axes[1].set_title("Family recovery")
    axes[1].set_ylabel("Median NMI")
    axes[2].set_title("Clan recovery")
    axes[2].set_ylabel("Median NMI")
    for axis in axes:
        axis.set_xlabel("Top fraction of available scores retained (%)")
        axis.legend()
    figure.suptitle(
        "Why equal score percentiles do not match BLAST and Biopython graphs"
    )
    figure.tight_layout()
    figure.savefig(output / "percentile_density_mismatch.png", dpi=200)
    plt.close(figure)


def plot_rule_frequency(
    results: pd.DataFrame,
    output: Path,
    *,
    metric: str,
    reference_clans: int,
    reference_families: int,
    reference_domains: int,
) -> None:
    """Plot one agreement distribution for all rules at one composition."""

    if metric not in {"nmi", "ami"}:
        raise ValueError("metric must be 'nmi' or 'ami'")

    subset = results.loc[
        (results["clans_per_collection"] == reference_clans)
        & (results["families_per_clan"] == reference_families)
        & (results["domains_per_family"] == reference_domains)
    ]
    if subset.empty:
        raise ValueError("the requested reference composition has no results")
    if not subset["analysis_eligible"].all():
        raise ValueError(
            "the requested reference composition contains repeated node sets; "
            "choose a primary-analysis-eligible cell"
        )

    settings = (
        subset[["threshold_rule", "threshold"]]
        .drop_duplicates()
        .sort_values(["threshold_rule", "threshold"])
    )
    bins = np.linspace(0, 1, 26)
    centres = (bins[:-1] + bins[1:]) / 2
    colours = plt.colormaps["tab20"](np.linspace(0, 1, len(settings)))
    linestyles = {"percentile": "-", "target_density": "--", "top_k": ":"}
    figure, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)
    for row, (reference, reference_title) in enumerate(REFERENCE_LABELS):
        for column, method in enumerate(("biopython", "blast")):
            axis = axes[row, column]
            panel = subset.loc[
                (subset["reference_label"] == reference) & (subset["method"] == method)
            ]
            for colour, setting in zip(colours, settings.itertuples(index=False)):
                values = panel.loc[
                    (panel["threshold_rule"] == setting.threshold_rule)
                    & (panel["threshold"] == setting.threshold),
                    metric,
                ]
                counts, _ = np.histogram(values, bins=bins)
                percentages = 100 * counts / len(values)
                axis.plot(
                    centres,
                    percentages,
                    color=colour,
                    linestyle=linestyles[setting.threshold_rule],
                    linewidth=1.8,
                    label=rule_label(setting.threshold_rule, setting.threshold),
                )
            axis.set_title(f"{method_label(method)} — {reference_title} labels")
            axis.set_xlabel(metric.upper())
            axis.set_ylabel(f"Collections in {metric.upper()} bin (%)")
            axis.set_xlim(0, 1)
            axis.grid(alpha=0.2)
    figure.suptitle(
        "All 12 graph rules at "
        f"{reference_clans} clans × {reference_families} families/clan × "
        f"{reference_domains} domains/family"
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="center left", bbox_to_anchor=(0.88, 0.5))
    figure.subplots_adjust(right=0.86)
    figure.tight_layout(rect=(0, 0, 0.86, 0.96))
    figure.savefig(
        output
        / (
            f"rule_{metric}_frequency_"
            f"c{reference_clans}_f{reference_families}_d{reference_domains}.png"
        ),
        dpi=200,
    )
    plt.close(figure)


def plot_composition_sweep(
    results: pd.DataFrame,
    uniqueness: pd.DataFrame,
    output: Path,
    *,
    method: str,
    rule: str,
    threshold: float,
) -> None:
    """Show one rule's NMI distribution in every composition cell."""

    subset = results.loc[
        (results["method"] == method)
        & (results["threshold_rule"] == rule)
        & (results["threshold"] == threshold)
        & (results["reference_label"] == "family_accession")
    ]
    if subset.empty:
        raise ValueError("the requested composition-sweep graph setting is absent")

    clans = sorted(subset["clans_per_collection"].unique())
    row_settings = [
        (family_count, domain_count)
        for family_count in sorted(subset["families_per_clan"].unique())
        for domain_count in sorted(subset["domains_per_family"].unique())
    ]
    status = uniqueness.set_index(COMPOSITION_COLUMNS)["analysis_eligible"]
    bins = np.linspace(0, 1, 21)
    centres = (bins[:-1] + bins[1:]) / 2
    figure, axes = plt.subplots(
        len(row_settings), len(clans), figsize=(13, 15), sharex=True, sharey=True
    )
    for row, (family_count, domain_count) in enumerate(row_settings):
        for column, clan_count in enumerate(clans):
            axis = axes[row, column]
            panel = subset.loc[
                (subset["clans_per_collection"] == clan_count)
                & (subset["families_per_clan"] == family_count)
                & (subset["domains_per_family"] == domain_count)
            ]
            eligible = bool(status.loc[(clan_count, family_count, domain_count)])
            counts, _ = np.histogram(panel["nmi"], bins=bins)
            percentages = 100 * counts / len(panel)
            colour = "#2878B5" if eligible else "#D07C2C"
            axis.fill_between(
                centres, percentages, step="mid", alpha=0.25, color=colour
            )
            axis.plot(centres, percentages, drawstyle="steps-mid", color=colour)
            axis.axvline(
                panel["nmi"].median(), color=colour, linestyle="--", linewidth=1
            )
            node_count = clan_count * family_count * domain_count
            axis.set_title(f"{clan_count} clans; n={node_count}")
            if column == 0:
                axis.set_ylabel(
                    f"{family_count} families/clan\n{domain_count} domains/family\nCollections (%)"
                )
            if row == len(row_settings) - 1:
                axis.set_xlabel("Family NMI")
            if not eligible:
                unique_sets = int(
                    uniqueness.loc[
                        (uniqueness["clans_per_collection"] == clan_count)
                        & (uniqueness["families_per_clan"] == family_count)
                        & (uniqueness["domains_per_family"] == domain_count),
                        "unique_node_sets",
                    ].iloc[0]
                )
                axis.set_facecolor("#FFF4E6")
                axis.text(
                    0.02,
                    0.94,
                    f"audit only: {unique_sets}/50 unique",
                    transform=axis.transAxes,
                    va="top",
                    fontsize=8,
                    color="#8A4B08",
                )
            axis.grid(alpha=0.15)
    figure.suptitle(
        f"Composition sweep — {method_label(method)}, {rule_label(rule, threshold)}, family recovery\n"
        "Blue cells enter distributional summaries; amber cells are retained for audit only"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    threshold_slug = f"{threshold:g}".replace(".", "p")
    figure.savefig(
        output / f"composition_sweep_{method}_{rule}_{threshold_slug}_family.png",
        dpi=200,
    )
    plt.close(figure)


def plot_nmi_ami_comparison(
    results: pd.DataFrame,
    output: Path,
    *,
    reference_clans: int,
    reference_families: int,
    reference_domains: int,
) -> None:
    """Compare raw and chance-adjusted partition agreement at one cell."""

    subset = results.loc[
        (results["analysis_eligible"])
        & (results["clans_per_collection"] == reference_clans)
        & (results["families_per_clan"] == reference_families)
        & (results["domains_per_family"] == reference_domains)
    ].copy()
    settings = (
        subset[["threshold_rule", "threshold"]]
        .drop_duplicates()
        .sort_values(["threshold_rule", "threshold"])
    )
    labels = [
        rule_label(row.threshold_rule, row.threshold) for row in settings.itertuples()
    ]
    x = np.arange(len(settings))
    figure, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True, sharey=True)
    for row, (reference, reference_title) in enumerate(REFERENCE_LABELS):
        for column, method in enumerate(("biopython", "blast")):
            axis = axes[row, column]
            panel = subset.loc[
                (subset["reference_label"] == reference) & (subset["method"] == method)
            ]
            medians = (
                panel.groupby(["threshold_rule", "threshold"])[["nmi", "ami"]]
                .median()
                .reindex(pd.MultiIndex.from_frame(settings))
            )
            axis.plot(x, medians["nmi"], marker="o", label="NMI")
            axis.plot(x, medians["ami"], marker="o", label="AMI")
            axis.fill_between(x, medians["ami"], medians["nmi"], alpha=0.12)
            axis.set_title(f"{method_label(method)} — {reference_title} labels")
            axis.set_ylabel("Median agreement")
            axis.set_ylim(-0.05, 1.02)
            axis.grid(alpha=0.2)
            axis.legend()
    for axis in axes[-1]:
        axis.set_xticks(x, labels, rotation=55, ha="right")
    figure.suptitle(
        "NMI versus chance-adjusted AMI at "
        f"{reference_clans} clans × {reference_families} families/clan × "
        f"{reference_domains} domains/family"
    )
    figure.tight_layout()
    figure.savefig(output / "nmi_vs_ami_by_rule.png", dpi=200)
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("full_results", type=Path)
    parser.add_argument("matched_density_results", type=Path)
    parser.add_argument("collection_manifest", type=Path)
    parser.add_argument("collection_membership", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("figure_directory", type=Path)
    parser.add_argument("--reference-clans", type=int, default=4)
    parser.add_argument("--reference-families-per-clan", type=int, default=2)
    parser.add_argument("--reference-domains-per-family", type=int, default=20)
    parser.add_argument(
        "--composition-method", choices=("biopython", "blast"), default="biopython"
    )
    parser.add_argument(
        "--composition-rule",
        choices=("top_k", "percentile", "target_density"),
        default="top_k",
    )
    parser.add_argument("--composition-threshold", type=float, default=5)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    full = pd.read_csv(arguments.full_results, sep="\t")
    matched = pd.read_csv(arguments.matched_density_results, sep="\t")
    collections = pd.read_csv(arguments.collection_manifest, sep="\t")
    membership = pd.read_csv(arguments.collection_membership, sep="\t")
    results = pd.concat([full, matched], ignore_index=True, sort=False)

    output = arguments.output_directory.resolve()
    figures = arguments.figure_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    fingerprints = collection_fingerprints(membership)
    collection_audit = collections.merge(
        fingerprints, on=["collection_id", "node_count"], validate="one_to_one"
    )
    uniqueness = collection_audit.groupby(
        ["clans_per_collection", "families_per_clan", "domains_per_family"],
        as_index=False,
    ).agg(
        declared_collections=("collection_id", "size"),
        unique_node_sets=("node_set_sha256", "nunique"),
    )
    uniqueness["duplicate_collection_count"] = (
        uniqueness["declared_collections"] - uniqueness["unique_node_sets"]
    )
    uniqueness["analysis_eligible"] = uniqueness["duplicate_collection_count"] == 0
    uniqueness["analysis_status"] = np.where(
        uniqueness["analysis_eligible"],
        "primary_independent_replicates",
        "audit_only_repeated_node_sets",
    )
    results, collection_status = add_analysis_eligibility(
        results, collections, uniqueness
    )
    collection_audit = collection_audit.merge(
        collection_status[["collection_id", "analysis_eligible"]],
        on="collection_id",
        how="left",
        validate="one_to_one",
    )
    primary_results = results.loc[results["analysis_eligible"]].copy()
    matched_target_densities = fully_matched_target_densities(results)

    rule_summary = quantile_summary(
        primary_results,
        ["reference_label", "method", "threshold_rule", "threshold"],
    )
    factor_summary = quantile_summary(
        primary_results,
        [
            "reference_label",
            "method",
            "threshold_rule",
            "threshold",
            "clans_per_collection",
            "families_per_clan",
            "domains_per_family",
            "node_count",
        ],
    )
    paired = paired_method_differences(primary_results)
    paired["comparison_status"] = paired.apply(
        lambda row: comparison_status(
            row["threshold_rule"], row["threshold"], matched_target_densities
        ),
        axis=1,
    )
    paired_summary = paired.groupby(
        ["reference_label", "threshold_rule", "threshold"], as_index=False
    ).agg(
        pairs=("collection_id", "size"),
        comparison_status=("comparison_status", "first"),
        median_nmi_difference=("nmi_blast_minus_biopython", "median"),
        median_absolute_nmi_difference=(
            "nmi_blast_minus_biopython",
            lambda values: values.abs().median(),
        ),
        lower_quartile_nmi_difference=(
            "nmi_blast_minus_biopython",
            lambda values: values.quantile(0.25),
        ),
        upper_quartile_nmi_difference=(
            "nmi_blast_minus_biopython",
            lambda values: values.quantile(0.75),
        ),
        median_density_difference=("density_blast_minus_biopython", "median"),
    )

    tables = {
        "collection_fingerprints.tsv": collection_audit,
        "collection_uniqueness.tsv": uniqueness,
        "rule_summary.tsv": rule_summary,
        "factor_summary.tsv": factor_summary,
        "paired_method_differences.tsv": paired,
        "paired_method_summary.tsv": paired_summary,
    }
    for filename, table in tables.items():
        table.to_csv(output / filename, sep="\t", index=False)

    plot_top_k_input_sensitivity(primary_results, figures)
    plot_matched_density_method_deltas(paired, figures, matched_target_densities)
    plot_percentile_mismatch(primary_results, figures)
    for metric in ("nmi", "ami"):
        plot_rule_frequency(
            primary_results,
            figures,
            metric=metric,
            reference_clans=arguments.reference_clans,
            reference_families=arguments.reference_families_per_clan,
            reference_domains=arguments.reference_domains_per_family,
        )
    plot_composition_sweep(
        results,
        uniqueness,
        figures,
        method=arguments.composition_method,
        rule=arguments.composition_rule,
        threshold=arguments.composition_threshold,
    )
    plot_nmi_ami_comparison(
        primary_results,
        figures,
        reference_clans=arguments.reference_clans,
        reference_families=arguments.reference_families_per_clan,
        reference_domains=arguments.reference_domains_per_family,
    )

    manifest_path = output / "summary_manifest.json"
    manifest = {
        "experiment": "full repeated Pfam graph experiment summary",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "graph_count": int(results["graph_id"].nunique()),
        "evaluation_row_count": len(results),
        "collection_count": int(results["collection_id"].nunique()),
        "primary_analysis_collection_count": int(
            primary_results["collection_id"].nunique()
        ),
        "excluded_composition_cell_count": int(
            (~uniqueness["analysis_eligible"]).sum()
        ),
        "analysis_rule": (
            "distributional summaries use only composition cells with 50/50 "
            "distinct node sets"
        ),
        "fully_matched_target_densities": matched_target_densities,
        "figure_configuration": {
            "reference_clans": arguments.reference_clans,
            "reference_families_per_clan": arguments.reference_families_per_clan,
            "reference_domains_per_family": arguments.reference_domains_per_family,
            "composition_method": arguments.composition_method,
            "composition_rule": arguments.composition_rule,
            "composition_threshold": arguments.composition_threshold,
        },
        "inputs": {
            str(path): sha256_file(path)
            for path in (
                arguments.full_results,
                arguments.matched_density_results,
                arguments.collection_manifest,
                arguments.collection_membership,
            )
        },
        "tables": {filename: sha256_file(output / filename) for filename in tables},
        "figures": {
            path.name: sha256_file(path) for path in sorted(figures.glob("*.png"))
        },
        "generator": "scripts/summarize_resampled_graph_experiment.py",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Summarized {manifest['graph_count']} graphs from "
        f"{manifest['collection_count']} collections"
    )
    print(f"Tables: {output}")
    print(f"Figures: {figures}")


if __name__ == "__main__":
    main()
