# Verification record (2026-10-22)

- All three notebooks executed cell by cell from top to bottom in an IPython session using a non-interactive Matplotlib backend. This did not use a Jupyter kernel/nbclient execution service.
- UVLF: 150 rendered observations, 55 primary observations, seven redshift panels and 14 cached theory curves; all notebook input, data-preservation, layout and exported-file assertions passed. Rendered PNG inspected.
- Accretion history: three central values with asymmetric intervals and three upper limits; notebook transform and export checks passed. Rendered PNG inspected.
- Posterior: all code cells completed with the portable snapshot. GetDist/corner emitted sparse-contour, smoothing or density warnings; full-chain autocorrelation was intentionally unavailable because the original chain was not installed inside this package. Completion is not a statistical validation.
- Snapshot arrays and log probabilities were checked for exact equality against emcee's original-chain output with discard=30000, thin=1000 and the recorded keep mask.
- MCMC script and exporter passed Python syntax parsing. The fresh/resume initial-state decision was checked; a new MPI sampling run and a new likelihood evaluation were not performed.
- Original supplied files were checked against their SHA-256 hashes and remain unchanged.
- Distribution notebooks contain no saved outputs or execution counts. No full HDF5 chain, .DS_Store, local provenance paths or Python caches are included in the upload list.
- The environment used existing UVLF dependencies plus getdist 1.7.7 and corner 2.3.0 in a temporary dependency directory, and the bundled pypdf package. This is not evidence of historical dependency versions or a clean-environment installation test.

Known scientific inconsistencies are listed in REPRODUCIBILITY.md and remain unresolved. The figures here are replay previews, not proof of current PRL formatting compliance. The posterior-derived figure was generated successfully but has not received a separate visual layout audit in this preparation.
