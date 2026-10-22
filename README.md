# PBH-UVLF

Analysis code and supporting data for a study of primordial black holes (PBHs) and galaxy ultraviolet luminosity functions over the seven redshift nodes z = 6, 7, 8, 9, 10, 11, and 12.5.

## What is included

| Path | Contents |
| --- | --- |
| `mcmc/` | Supplied 16-parameter model and MPI/emcee sampler; UVLF likelihood input and dust-calibration table |
| `plot/UVLF_full_WithPBH_noRobertson_20261016.ipynb` | Seven-panel UVLF plot using the supplied cached theoretical curves |
| `plot/prl_free_sf_accretion_rate_v22.ipynb` | Plot of six supplied Eddington-ratio constraints and corresponding accretion rates |
| `plot/statistics and plot_F_216_consistent_MCSamples_shortest68.ipynb` | Supplied posterior-analysis notebook, adapted to read the portable sample file |
| `plot/*.txt` | Literature observation ledgers, cached theory curves, and plotting inputs/exports |
| `posterior/posterior_samples.npz` | 6512 equally weighted samples selected using the supplied notebook's processing rules |
| `posterior/metadata.json` | Chain checksum, processing settings, excluded walker indices and array dimensions |
| `posterior/parameters.csv` | Column order and definitions for the posterior arrays |
| `scripts/export_posterior.py` | Recreate the portable posterior file from the full HDF5 chain |
| `figures/` | PDF previews regenerated from the supplied notebooks |
| `docs/REPRODUCIBILITY.md` | Verification scope, known inconsistencies and changes made for distribution |
| `docs/SOURCE_MANIFEST.json` | Provenance and original checksums of supplied or recovered inputs |

This release preserves the supplied scientific calculations. It is not a new MCMC run or a validation of the paper's statistical claims. Review `docs/REPRODUCIBILITY.md` before using the posterior intervals.

## Installation

Use Python 3.12 or a compatible newer version in a dedicated environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The version list records the environment used for the release checks; it does not reconstruct the historical chain-generation environment. The full fresh-install procedure has not been tested.

## Recreate the plots

```bash
cd plot
jupyter lab
```

Open either plotting notebook and run its cells from top to bottom. Keep the notebook working directory at `plot/`; filenames in these notebooks are relative to that directory. The UVLF notebook reads cached model curves and does not perform a new fit. The accretion-rate notebook plots supplied numerical constraints and does not infer them from the posterior snapshot. Times New Roman is requested; systems without that font may render a fallback.

The posterior-analysis notebook reads `../posterior/posterior_samples.npz`. Original processing: discard 30000 steps, reject walkers whose post-burn-in fraction of moving steps is below 0.02, then use emcee thinning by 1000. There are 74 retained draws per kept walker, 88 kept walkers and 16 parameters. Full-resolution traces and autocorrelation diagnostics require the original chain; thinned traces must not be interpreted as full-chain convergence evidence.

## Full MCMC sampling (optional, expensive)

Install an MPI implementation (e.g. Open MPI) and the optional dependencies:

```bash
python -m pip install -r requirements-mpi.txt
mpiexec -n 4 python mcmc/emcee_UVLF_P18_a16_merged_z20.py
```

The script uses 112 walkers and 60000 steps per invocation, with the original prior and model settings. A fresh HDF5 backend starts from the supplied initialization ranges; an existing backend resumes. No existing chain is reset. The output is `mcmc/emcee_UVLF_P18_a16_merged_z20.h5`. The script does not fix a random seed, so a fresh chain will not reproduce historical samples byte for byte. This expensive sampling command was not executed for the release.

The full original chain is about 1.9 GB and is not included in this GitHub repository. Its SHA-256 and selection details are recorded in `posterior/metadata.json`; this repository provides the processed sample snapshot, not full-chain diagnostic data. No public download location for the full chain is asserted.

## Data provenance and units

The fitting input `mcmc/UVLF_PBH.txt` has six columns: redshift; M_UV (AB mag); magnitude-bin width (mag); Phi_UV (Mpc^-3 mag^-1); upper uncertainty; lower uncertainty. The supplied likelihood applies a minimum 20% uncertainty and its own dust/cosmology transformations; keep this input table unmodified. Its header identifies the HST source (arXiv:2102.07775), but the high-redshift block needs source confirmation; see the reproducibility notes.

`mcmc/Beta_parameters.txt` records the dust UV-slope calibration and identifies arXiv:1306.2950, Table 3. `plot/UVLF_literature_data.txt` and the exported `*_data.txt` provide row-level references, source tables, redshift intervals, uncertainty definitions, limits and display flags. Comparison-only points are not additional likelihood inputs. Published observational cosmologies and selection functions have not all been homogenized.

Theory tables record M_UV in AB mag and number densities in Mpc^-3 mag^-1, with `phi_model` and `phi_observed` distinguished in their headers. Posterior logarithms are base 10; `log_prob` is the emcee log posterior, not a separately recomputed log likelihood. See `posterior/parameters.csv` for parameter order. Accretion-rate conversion in the supplied figure is mdot = 9.73 lambda_E.

## Citation and availability

Repository: https://github.com/hudiao1/PBH-UVLF

A manuscript citation or archive DOI has not been added because no verified publication identifier or DOI was supplied. Cite the original observational publications in addition to this repository when reusing their measurements. A tagged release and a persistent archive can be added when the authors freeze the paper version.

## Rights and licensing

No license has been selected for this repository. Public access does not grant a general reuse or redistribution license. Rights in author-created code remain with their respective rights holders; contact them for reuse permissions. Third-party data and software retain their original rights and conditions. Source citations document provenance and are not substitutes for permission. Do not assume the compiled literature tables have been relicensed.
