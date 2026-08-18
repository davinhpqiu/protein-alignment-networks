"""Run pretrained DEDAL inference in its isolated TensorFlow environment."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import tensorflow_hub as hub
from dedal import infer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--model", required=True)
    arguments = parser.parse_args()

    model = hub.load(arguments.model)
    with arguments.input.open(encoding="utf-8") as source, arguments.output.open(
        "w", encoding="utf-8"
    ) as destination:
        for line in source:
            pair = json.loads(line)
            result = align_pair(
                model,
                pair["protein_a"],
                pair["protein_b"],
                pair["sequence_a"],
                pair["sequence_b"],
            )
            destination.write(json.dumps(result) + "\n")
            destination.flush()


def align_pair(model, protein_a: str, protein_b: str, sequence_a: str, sequence_b: str):
    inputs = infer.preprocess(sequence_a, sequence_b)
    output = model(inputs)
    expanded = infer.expand(
        [output["sw_scores"], output["paths"], output["sw_params"]]
    )
    processed = infer.postprocess(expanded, len(sequence_a), len(sequence_b))
    alignment = infer.Alignment(sequence_a, sequence_b, *processed)
    logit = float(output["homology_logits"].numpy()[0])
    probability = 1.0 / (1.0 + math.exp(-logit)) if logit >= 0 else (
        math.exp(logit) / (1.0 + math.exp(logit))
    )
    return {
        "protein_a": protein_a,
        "protein_b": protein_b,
        "sw_score": float(output["sw_scores"].numpy()[0]),
        "homology_logit": logit,
        "homology_probability": probability,
        "aligned_a": alignment.left_match,
        "alignment_symbols": alignment.matches,
        "aligned_b": alignment.right_match,
        "identity_count": int(alignment.identity),
        "similarity_count": int(alignment.similarity),
        "gap_count": int(alignment.gaps),
        "alignment_length": len(alignment),
        "a_start": int(alignment.start[0]),
        "a_end": int(alignment.end[0]),
        "b_start": int(alignment.start[1]),
        "b_end": int(alignment.end[1]),
    }


if __name__ == "__main__":
    main()
