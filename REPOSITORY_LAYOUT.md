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

## Vendored dependencies

`booksim2/` is vendored in this repository from BookSim 2.0 upstream commit
`28f43299f1706a3160ffac721ca461d74eb6e618`. Its CCDG and WSE extensions are
tracked directly by the main repository so one clone contains the complete
simulator implementation. Build products and simulation outputs remain ignored.

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
