from scripts.prepare_pfam_clan_panel import aligned_distance, select_diverse


def candidate(node_id, aligned_sequence):
    sequence = aligned_sequence.replace("-", "")
    return {
        "node_id": node_id,
        "domain_length": len(sequence),
        "aligned_sequence": aligned_sequence,
    }


def test_aligned_distance_ignores_gap_only_positions():
    assert aligned_distance("AC-D", "ACED") == 0.0
    assert aligned_distance("AAAA", "CCCC") == 1.0


def test_diverse_selection_is_deterministic_and_retains_extremes():
    candidates = [
        candidate("a", "AAAA"),
        candidate("b", "AAAC"),
        candidate("c", "CCCC"),
    ]
    first = select_diverse(candidates, 2)
    second = select_diverse(candidates, 2)
    assert [row["node_id"] for row in first] == [row["node_id"] for row in second]
    assert {row["node_id"] for row in first} == {"a", "c"}
