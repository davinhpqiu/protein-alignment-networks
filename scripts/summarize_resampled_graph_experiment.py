"""Summarize the full repeated Pfam graph experiment and create core figures."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    required = {(metric, method) for metric in metrics for method in ("biopython", "blast")}
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


def plot_top_k_input_sensitivity(results: pd.DataFrame, output: Path) -> None:
    subset = results.loc[
        (results["method"] == "biopython")
        & (results["threshold_rule"] == "top_k")
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


def plot_matched_density_method_deltas(paired: pd.DataFrame, output: Path) -> None:
    subset = paired.loc[paired["threshold_rule"] == "target_density"]
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
        axis.boxplot(values, tick_labels=[f"{100*x:g}%" for x in thresholds], showfliers=False)
        axis.axhline(0, color="black", linewidth=1)
        axis.set_title(title)
        axis.set_xlabel("Matched target edge density")
        axis.set_ylabel("BLAST NMI − Biopython NMI")
    figure.suptitle("Paired method differences after matching graph density")
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
    figure.suptitle("Why equal score percentiles do not match BLAST and Biopython graphs")
    figure.tight_layout()
    figure.savefig(output / "percentile_density_mismatch.png", dpi=200)
    plt.close(figure)


def plot_fixed_condition_nmi_distributions(results: pd.DataFrame, output: Path) -> None:
    subset = results.loc[
        (results["threshold_rule"] == "target_density")
        & (results["clans_per_collection"] == 4)
        & (results["families_per_clan"] == 2)
        & (results["domains_per_family"] == 20)
    ]
    bins = np.linspace(0, 1, 31)
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    for row, reference in enumerate(("family_accession", "clan_accession")):
        for column, method in enumerate(("biopython", "blast")):
            axis = axes[row, column]
            panel = subset.loc[
                (subset["reference_label"] == reference)
                & (subset["method"] == method)
            ]
            for threshold, group in panel.groupby("threshold"):
                axis.hist(
                    group["nmi"],
                    bins=bins,
                    density=True,
                    histtype="step",
                    linewidth=2,
                    label=f"density={100*threshold:g}%",
                )
            axis.set_title(f"{method.title()} — {reference.replace('_accession', '')}")
            axis.set_xlabel("NMI")
            axis.set_ylabel("Empirical density")
            axis.set_xlim(0, 1)
            axis.legend()
    figure.suptitle(
        "Matched-density NMI distributions: 4 clans × 2 families × 20 domains"
    )
    figure.tight_layout()
    figure.savefig(output / "matched_density_nmi_distributions.png", dpi=200)
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("full_results", type=Path)
    parser.add_argument("matched_density_results", type=Path)
    parser.add_argument("collection_manifest", type=Path)
    parser.add_argument("collection_membership", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("figure_directory", type=Path)
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
    uniqueness = (
        collection_audit.groupby(
            ["clans_per_collection", "families_per_clan", "domains_per_family"],
            as_index=False,
        )
        .agg(
            declared_collections=("collection_id", "size"),
            unique_node_sets=("node_set_sha256", "nunique"),
        )
    )
    uniqueness["duplicate_collection_count"] = (
        uniqueness["declared_collections"] - uniqueness["unique_node_sets"]
    )

    rule_summary = quantile_summary(
        results,
        ["reference_label", "method", "threshold_rule", "threshold"],
    )
    factor_summary = quantile_summary(
        results,
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
    paired = paired_method_differences(results)
    paired_summary = (
        paired.groupby(
            ["reference_label", "threshold_rule", "threshold"], as_index=False
        )
        .agg(
            pairs=("collection_id", "size"),
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

    plot_top_k_input_sensitivity(results, figures)
    plot_matched_density_method_deltas(paired, figures)
    plot_percentile_mismatch(results, figures)
    plot_fixed_condition_nmi_distributions(results, figures)

    manifest_path = output / "summary_manifest.json"
    manifest = {
        "experiment": "full repeated Pfam graph experiment summary",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "graph_count": int(results["graph_id"].nunique()),
        "evaluation_row_count": len(results),
        "collection_count": int(results["collection_id"].nunique()),
        "inputs": {
            str(path): sha256_file(path)
            for path in (
                arguments.full_results,
                arguments.matched_density_results,
                arguments.collection_manifest,
                arguments.collection_membership,
            )
        },
        "tables": {
            filename: sha256_file(output / filename) for filename in tables
        },
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
