# Reproducibility notes

## Distribution changes

- Original files in the author's `data_availability` directory were not edited.
- Notebook outputs and execution counters were cleared; code-cell order and scientific values were preserved.
- Missing UVLF inputs were recovered from the existing UVLF figure workspace. Their original checksums are in `SOURCE_MANIFEST.json`.
- The posterior notebook now reads a portable NPZ snapshot rather than requiring a 1.9 GB HDF5 file. The export uses precisely the supplied notebook's discard=30000, thin=1000 and movement<0.02 walker-selection rule. Acceptance rates are exported from the original chain. Full-chain autocorrelation is only attempted if the original HDF5 file is supplied under `mcmc/`.
- The MCMC script resolves its two input tables and output backend relative to its own directory.
- A sampling startup error was repaired: the supplied script computed initial positions but always passed None to `run_mcmc`. The release passes those positions for a fresh backend and None only when resuming an existing chain. Model functions, prior bounds, initialization ranges, walker count and step count were preserved.
- No new likelihood evaluation, fit, MCMC sampling, external-constraint reconstruction or chain-convergence claim is made by this distribution work.

## Scientific issues retained for author review

1. The original posterior notebook's section labelled `16%-84%` actually calls `np.percentile(..., [5,32,95])`. These values were not changed; this cell must not be described as a median and central 68% interval.
2. Several KDE cells use lambda_E hard bounds (-15,5), whereas the sampler and the Eddington triangle-plot cell use (-15,2). KDE smoothing choices also differ between cells. The authors must choose and verify a common reporting convention before treating all cells as one consistent table-generation procedure.
3. A mass/abundance plotting cell assigns q16,q50,q84 using [0.16,0.50,0.95]. Its variable names do not make this a 16-84% interval.
4. The accretion-rate plot contains author-supplied summary values directly. These values were not regenerated from this posterior snapshot or verified against manuscript Table I during this release preparation.
5. The MCMC input header explicitly identifies arXiv:2102.07775, but the file also contains z=9-12.5 points without a distinct block-level source citation. The manuscript identifies Donnan et al. (2024) for high-redshift inputs; exact table/version/row correspondence still requires author verification.
6. Cached UVLF theory files record historical source hashes. The supplied sampler file was not proven to be that exact historical generating version. Cached curve replay is distinct from regenerating the curves from this sampler.
7. The supplied plotting notebook describes a readability variant with 1 mm markers and a 5 pt PDF text floor; reproducing it is not certification of PRL figure-format compliance.
8. Posterior samples are selected using the supplied notebook's walker exclusion rule. This release does not establish the validity of that exclusion or MCMC convergence. The original chain has 104666 stored steps, while the supplied script uses 60000 steps per invocation; the continuation history has not been reconstructed.

## Full-chain data

The original HDF5 chain is excluded from ordinary Git tracking. The processed NPZ file is suitable for replaying the preserved posterior-analysis cells, subject to the issues above. It cannot reconstruct rejected walkers, all intermediate steps or full-chain autocorrelation. Ask the authors about access to the original chain; this release does not commit them to a specific future download route.

The compact snapshot includes chains[draw,walker,parameter], log_prob[draw,walker], original-walker keep mask, excluded walker indices, full post-burn-in movement fractions, accepted counts, iteration count, discard and thinning. Sample weights are equal. The exporter records a full-file SHA-256 of the source chain.

## Verification

See `VERIFICATION.md` for the checks actually run. Static syntax checks and plotting replay do not validate the scientific assumptions or reported significance.
