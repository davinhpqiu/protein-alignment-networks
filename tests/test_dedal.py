import math

import pytest

from protein_alignment_networks.dedal import dedal_all_vs_all, homology_probability


def test_homology_probability_is_stable():
    assert homology_probability(0.0) == 0.5
    assert homology_probability(1000.0) == pytest.approx(1.0)
    assert homology_probability(-1000.0) == pytest.approx(0.0)
    assert homology_probability(math.log(3)) == pytest.approx(0.75)


def test_dedal_rejects_sequences_that_would_be_cropped():
    with pytest.raises(ValueError, match="at most 511"):
        dedal_all_vs_all({"short": "MKT", "long": "A" * 512})


def test_dedal_single_sequence_needs_no_model():
    result = dedal_all_vs_all({"only": "MKT"})
    assert result.empty
