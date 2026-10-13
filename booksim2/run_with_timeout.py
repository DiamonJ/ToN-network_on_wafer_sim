#!/usr/bin/env python3
"""Run a command with a portable wall-clock timeout.

GNU ``timeout`` is not installed by default on macOS.  The BookSim runners
use this helper so a timeout is exit status 124 instead of the misleading 127
(``command not found``).
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys


def main() -> int:
    if len(sys.argv) < 3:
        print(f"usage: {sys.argv[0]} <seconds> <command> [args ...]", file=sys.stderr)
        return 2
    try:
        seconds = float(sys.argv[1])
    except ValueError:
        print(f"invalid timeout: {sys.argv[1]}", file=sys.stderr)
        return 2
    if seconds <= 0:
        print("timeout must be positive", file=sys.stderr)
        return 2

    proc = subprocess.Popen(sys.argv[2:], start_new_session=True)
    try:
        return proc.wait(timeout=seconds)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
        return 124


if __name__ == "__main__":
    raise SystemExit(main())
