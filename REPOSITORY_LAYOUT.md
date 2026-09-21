# Repository layout and artifact policy

This repository keeps project-owned source, configuration, tests, small input
cases, and compact result summaries in Git. Large or reproducible build and run
artifacts stay outside Git history.

## Tracked in the main repository

- Top-level Python and shell pipeline tools.
- `tests/`, `docs/`, `simgrid_traces/`, and analysis source under
  `traffic_pattern/`.
- Small LAMMPS inputs and potential files under `cases/`.
- Project-owned LAMMPS instrumentation files under `lammps-src/src/`. The full
  upstream LAMMPS source tree is a local dependency and is not vendored here.
- Compact experiment summaries (`csv`, report `json`, and selected figures).

## Separate dependency repositories

`booksim2/` is an independent Git checkout based on BookSim 2.0. Its CCDG and
WSE extensions are committed in that repository so upstream history remains
available. The main repository intentionally ignores the nested checkout.

`lammps-src/` and `sst-dumpi/` are upstream source trees used to build the
pipeline. Keep a known-compatible checkout locally, then apply or copy the
project-owned instrumentation files tracked by this repository.

## Generated locally

The following paths are intentionally ignored because they are reproducible or
machine-specific:

- `install/` and native build products.
- `runs/`, DUMPI traces, logs, and full CCDG/WSE replay payloads.
- Generated `*.program.json`, `*.replay.ccdg`, BookSim configs, and stats.
- Reproduction archives such as `lammps_trace_repro_pkg.tar.gz`.

Publish large, non-reproducible artifacts in object storage or a release and
record their URL, checksum, generating commit, command, and configuration in a
small manifest. Git LFS should be reserved for essential binary inputs that
cannot be regenerated.
