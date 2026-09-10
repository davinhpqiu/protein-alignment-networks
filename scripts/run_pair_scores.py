"""Run the configurable all-pairs scoring pipeline for any prepared dataset.

This dataset-neutral entry point delegates to the shared runner historically
named ``run_pf00042_correlation.py``.  It writes method-specific outputs, one
canonical merged pair table, correlations, and a manifest containing parameters,
versions, runtimes, and checksums.  The notebooks load those artifacts rather
than maintaining copied numerical results.
"""

from run_pf00042_correlation import main

if __name__ == "__main__":
    main()
