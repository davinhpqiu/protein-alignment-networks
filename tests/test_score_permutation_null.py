import pandas as pd

from scripts.run_score_permutation_null import permute_finite_scores, stable_seed


def test_score_permutation_is_reproducible_and_preserves_missingness():
    pairs = pd.DataFrame(
        {
            "protein_a": ["a", "a", "b", "b"],
            "protein_b": ["b", "c", "c", "d"],
            "score": [1.0, 2.0, float("nan"), 4.0],
        }
    )
    first = permute_finite_scores(pairs, "score", seed=7)
    second = permute_finite_scores(pairs, "score", seed=7)

    pd.testing.assert_frame_equal(first, second)
    assert first["score"].isna().tolist() == pairs["score"].isna().tolist()
    assert sorted(first["score"].dropna()) == sorted(pairs["score"].dropna())
    assert first[["protein_a", "protein_b"]].equals(
        pairs[["protein_a", "protein_b"]]
    )


def test_stable_seed_changes_with_condition():
    assert stable_seed(1, "a", 5) == stable_seed(1, "a", 5)
    assert stable_seed(1, "a", 5) != stable_seed(1, "a", 6)
