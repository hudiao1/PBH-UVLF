"""Export the existing notebook's posterior selection without rerunning MCMC."""
import argparse
import hashlib
import json
from pathlib import Path
import emcee
import h5py
import numpy as np


def export(chain_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    reader = emcee.backends.HDFBackend(str(chain_path), read_only=True)
    discard, thin = 15000, 1000
    if reader.iteration <= discard + thin:
        raise ValueError('Chain is too short for the original notebook selection.')
    with h5py.File(chain_path, 'r') as f:
        dataset = f['mcmc/chain']
        counts = np.zeros(dataset.shape[1], dtype=np.int64)
        previous = dataset[discard]
        for start in range(discard + 1, reader.iteration, 4096):
            block = dataset[start:min(start + 4096, reader.iteration)]
            counts += np.any(np.diff(np.concatenate((previous[None], block)), axis=0) != 0, axis=2).sum(axis=0)
            previous = block[-1]
    move_fraction = counts / (reader.iteration - discard - 1)
    bad_walkers = np.flatnonzero(move_fraction < 0.02)
    keep = np.ones(reader.shape[0], dtype=bool)
    keep[bad_walkers] = False
    if not keep.any():
        raise ValueError('The original movement filter removed every walker.')
    chains = reader.get_chain(discard=discard, thin=thin)[:, keep, :]
    log_prob = reader.get_log_prob(discard=discard, thin=thin)[:, keep]
    assert chains.shape[:2] == log_prob.shape
    assert np.isfinite(chains).all() and np.isfinite(log_prob).all()
    np.savez_compressed(output_dir / 'posterior_samples.npz', chains=chains,
        log_prob=log_prob, keep=keep, bad_walkers=bad_walkers,
        move_fraction=move_fraction, accepted=reader.accepted,
        iteration=reader.iteration, discard=discard, thin=thin)
    digest = hashlib.sha256()
    with chain_path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(block)
    metadata = dict(source_file=chain_path.name, source_sha256=digest.hexdigest(),
        iteration=int(reader.iteration), original_walkers=int(reader.shape[0]),
        discard=discard, thin=thin, movement_threshold=0.02,
        excluded_walker_indices=bad_walkers.tolist(), retained_shape=list(chains.shape),
        flattened_samples=int(np.prod(chains.shape[:2])), weights='equal',
        log_prob_definition='emcee log posterior; not relabelled as log likelihood',
        selection='Identical to the supplied notebook: full post-burn-in walker movement filter, then emcee thinning and flattening.',
        limitation='This snapshot cannot establish convergence or reproduce full-resolution traces/autocorrelation. No new MCMC run was performed.')
    (output_dir / 'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('chain', type=Path)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'posterior')
    args = parser.parse_args()
    export(args.chain, args.output)
