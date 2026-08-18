"""Local adapter for the pretrained DEDAL research model."""

from __future__ import annotations

import json
import math
import os
import subprocess
import tempfile
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from .io import validate_protein_sequence

DEDAL_COLUMNS = [
    "protein_a",
    "protein_b",
    "sw_score",
    "homology_logit",
    "homology_probability",
    "aligned_a",
    "alignment_symbols",
    "aligned_b",
    "identity_count",
    "similarity_count",
    "gap_count",
    "alignment_length",
    "a_start",
    "a_end",
    "b_start",
    "b_end",
]


def dedal_all_vs_all(
    sequences: Mapping[str, str],
    *,
    output_path: str | Path | None = None,
    python_executable: str | Path | None = None,
    model: str | Path | None = None,
    project_root: str | Path | None = None,
) -> pd.DataFrame:
    """Align every unique pair with the locally cached pretrained DEDAL model.

    DEDAL is executed in its isolated environment so TensorFlow does not become
    a dependency of the main package. Sequences are limited to 511 residues:
    the published model input has 512 positions including its EOS token.
    """

    if not sequences:
        raise ValueError("at least one sequence is required")
    clean = {
        identifier: validate_protein_sequence(sequence, allow_empty=False)
        for identifier, sequence in sequences.items()
    }
    too_long = [identifier for identifier, sequence in clean.items() if len(sequence) > 511]
    if too_long:
        raise ValueError(
            "DEDAL supports at most 511 residues plus EOS; too long: "
            + ", ".join(too_long)
        )

    root = Path(project_root) if project_root else Path(__file__).resolve().parents[2]
    python = Path(python_executable) if python_executable else root / ".venv-dedal/bin/python"
    worker = root / "scripts/dedal_inference.py"
    source = root / "external/google-research"
    cache = root / "models/dedal"
    model_handle = str(model or "https://tfhub.dev/google/dedal/3")

    for required, description in [
        (python, "DEDAL Python environment"),
        (worker, "DEDAL inference worker"),
        (source / "dedal", "Google Research DEDAL source"),
    ]:
        if not required.exists():
            raise RuntimeError(f"{description} not found at {required}; see docs/DEDAL.md")

    names = list(clean)
    pairs = [
        {
            "protein_a": names[i],
            "protein_b": names[j],
            "sequence_a": clean[names[i]],
            "sequence_b": clean[names[j]],
        }
        for i in range(len(names))
        for j in range(i + 1, len(names))
    ]
    if not pairs:
        return pd.DataFrame(columns=DEDAL_COLUMNS)

    with tempfile.TemporaryDirectory(prefix="protein_alignment_dedal_") as temporary:
        temporary_path = Path(temporary)
        input_path = temporary_path / "pairs.jsonl"
        result_path = temporary_path / "results.jsonl"
        input_path.write_text(
            "".join(json.dumps(pair) + "\n" for pair in pairs),
            encoding="utf-8",
        )
        environment = os.environ.copy()
        environment["TFHUB_CACHE_DIR"] = str(cache)
        environment["PYTHONPATH"] = os.pathsep.join(
            [str(source), environment.get("PYTHONPATH", "")]
        ).rstrip(os.pathsep)
        _run_worker(
            [
                str(python),
                str(worker),
                "--input",
                str(input_path),
                "--output",
                str(result_path),
                "--model",
                model_handle,
            ],
            environment,
        )
        rows = [json.loads(line) for line in result_path.read_text().splitlines() if line]

    results = pd.DataFrame(rows, columns=DEDAL_COLUMNS)
    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(destination, sep="\t", index=False)
    return results


def homology_probability(logit: float) -> float:
    """Convert a DEDAL homology logit to its logistic probability."""

    if logit >= 0:
        return 1.0 / (1.0 + math.exp(-logit))
    exponential = math.exp(logit)
    return exponential / (1.0 + exponential)


def _run_worker(command: list[str], environment: dict[str, str]) -> None:
    try:
        subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
    except FileNotFoundError as error:
        raise RuntimeError(f"DEDAL executable not found: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip() or "no diagnostic output"
        raise RuntimeError(f"DEDAL inference failed: {detail}") from error
