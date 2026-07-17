#!/usr/bin/env python3
"""Parse AMPS reference-data dump files written by the instrumented
scale_atmos_phy_mp_amps.F90 (branch cloudlab_port).

Binary layout v1: sequence of records; ints int32, reals float64, byte order
auto-detected per file (little-endian by default; some cluster builds force
big-endian unformatted/stream I/O, e.g. gnu Makedefs' -fconvert=big-endian
or intel's -convert big_endian — this applies to stream writes too), arrays
prefixed by one int32 size per rank, Fortran (column-major) order. Record
types: micro (magic 1095586131) and sed (magic 1095586132) — field order
defined in the M0 plan Tasks 2-3 and mirrored here.

Every rank writes its own files, named amps_dump_r{RRRRRR}_t{TTT}.bin (rank,
thread); every rank dumps the same LOCAL i/j box, so npz keys must be
rank-qualified (see `_key`/`aggregate`) or cross-rank aggregation silently
collides and drops data.

Usage: python amps_dump_reader.py DUMP_DIR -o out.npz
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

MAGIC_MICRO = 1095586131
MAGIC_SED = 1095586132

MICRO_HEADER = ("phase", "TIME_AMPS", "i", "j", "isect", "nmic",
                "npr", "nbr", "ncr", "npi", "nbi", "nci", "npa", "nba", "nca",
                "mxnbin", "istrt", "jseed", "ifrst", "isect_seed", "nextn")
MICRO_R1 = ("qcvm", "v3v", "qvvm", "moist_denvm", "ptotvm", "tvm", "wbvm",
            "trpv_thil", "trpv_qtp")
MICRO_R4 = ("qrpvm", "qipvm", "qapvm")
MICRO_POST_EXTRA = (("dmtendlm", 3), ("dcontendlm", 3), ("dbintendlm", 4))

SED_HEADER = ("phase", "TIME_AMPS", "i", "j", "isect", "isn",
              "iadvv", "np", "nb", "nc", "k1", "k2", "k1m", "k2m")
SED_R1 = ("q_this", "q_other", "qcv", "qtp", "moist_denv", "thetav", "qvv",
          "tv", "dens_col", "momz_col", "u_col", "v_col", "cz_col", "fz_col",
          "dzzmv", "dzvmv")
SED_R1_TAIL = ("den_t", "momz_t", "rhou_t", "rhov_t", "rhoe_t")


class _Cursor:
    def __init__(self, buf: bytes, bo: str = "<"):
        self.buf = buf
        self.pos = 0
        self.bo = bo  # "<" little-endian (default) or ">" big-endian; one
        # value per file, auto-detected from the first record's magic in
        # read_dump_file (all records in a file share the build's byte order).

    def eof(self) -> bool:
        return self.pos >= len(self.buf)

    def i4(self, n: int = 1):
        out = np.frombuffer(self.buf, dtype=f"{self.bo}i4", count=n, offset=self.pos)
        self.pos += 4 * n
        return int(out[0]) if n == 1 else out.copy()

    def f8(self, n: int = 1):
        out = np.frombuffer(self.buf, dtype=f"{self.bo}f8", count=n, offset=self.pos)
        self.pos += 8 * n
        return float(out[0]) if n == 1 else out.copy()

    def i1_arr(self) -> np.ndarray:
        n = self.i4()
        return self.i4(n) if n > 0 else np.empty(0, dtype=f"{self.bo}i4")

    def rn_arr(self, ndims: int) -> np.ndarray:
        shape = tuple(self.i4() for _ in range(ndims))
        n = int(np.prod(shape))
        flat = self.f8(n) if n > 0 else np.empty(0)
        return np.asarray(flat).reshape(shape, order="F")


def _read_micro(c: _Cursor) -> dict:
    rec: dict = {"kind": "micro"}
    version = c.i4()
    assert version == 1, f"unsupported micro record version {version}"
    for name in MICRO_HEADER:
        rec[name] = c.i4()
    rec["dt"] = c.f8()
    rec["kmicvm"] = c.i1_arr()
    for name in MICRO_R1:
        rec[name] = c.rn_arr(1)
    for name in MICRO_R4:
        rec[name] = c.rn_arr(4)
    if rec["phase"] == 2:
        for name, ndims in MICRO_POST_EXTRA:
            rec[name] = c.rn_arr(ndims)
    return rec


def _read_sed(c: _Cursor) -> dict:
    rec: dict = {"kind": "sed"}
    version = c.i4()
    assert version == 1, f"unsupported sed record version {version}"
    for name in SED_HEADER:
        rec[name] = c.i4()
    rec["dt"] = c.f8()
    # k1b/k2b are written flattened (one value per (bin, class) column);
    # restore the Fortran (nb, nc) shape using the header's own nb/nc.
    rec["k1b"] = c.i1_arr().reshape((rec["nb"], rec["nc"]), order="F")
    rec["k2b"] = c.i1_arr().reshape((rec["nb"], rec["nc"]), order="F")
    rec["qpv"] = c.rn_arr(4)
    for name in SED_R1:
        rec[name] = c.rn_arr(1)
    rec["mmass"] = c.rn_arr(3)
    for name in SED_R1_TAIL:
        rec[name] = c.rn_arr(1)
    rec["sflx"] = c.f8()
    return rec


def _detect_byte_order(buf: bytes, path) -> str:
    """Return '<' or '>' for the file's byte order, detected from whichever
    interpretation of the first 4 bytes yields a known record magic. SCALE
    cluster builds commonly force big-endian unformatted/stream I/O
    (gnu -fconvert=big-endian, intel -convert big_endian), so this cannot be
    hardcoded. Every record in a file shares the build's byte order, so this
    is only done once, on the first magic."""
    raw = int(np.frombuffer(buf, dtype="<i4", count=1, offset=0)[0])
    if raw in (MAGIC_MICRO, MAGIC_SED):
        return "<"
    swapped = int(np.frombuffer(buf, dtype=">i4", count=1, offset=0)[0])
    if swapped in (MAGIC_MICRO, MAGIC_SED):
        return ">"
    raise ValueError(f"{path}: bad magic {raw} at byte 0")


def read_dump_file(path: str | Path) -> list[dict]:
    buf = Path(path).read_bytes()
    bo = "<" if len(buf) < 4 else _detect_byte_order(buf, path)
    c = _Cursor(buf, bo)
    records = []
    while not c.eof():
        magic = c.i4()
        if magic == MAGIC_MICRO:
            records.append(_read_micro(c))
        elif magic == MAGIC_SED:
            records.append(_read_sed(c))
        else:
            raise ValueError(f"{path}: bad magic {magic} at byte {c.pos - 4}")
    return records


_FNAME_RE = re.compile(r"amps_dump_r(\d+)_t(\d+)\.bin$")


def _parse_fname(path: Path) -> tuple[int, int]:
    """Extract (rank, thread) from an amps_dump_r{RRRRRR}_t{TTT}.bin name."""
    m = _FNAME_RE.search(path.name)
    if not m:
        raise ValueError(f"{path}: filename doesn't match amps_dump_r{{rank}}_t{{thread}}.bin")
    return int(m.group(1)), int(m.group(2))


def _key(rec: dict, rank: int) -> str:
    # Thread is deliberately NOT part of the key: within one rank, a given
    # (i, j) column is processed by exactly one thread, so (rank, TIME_AMPS,
    # i, j) is already unique. Rank IS required: every rank dumps the same
    # LOCAL i/j box, so without it, keys collide across ranks.
    if rec["kind"] == "micro":
        ph = "pre" if rec["phase"] == 1 else "post"
        return f"micro_r{rank}_t{rec['TIME_AMPS']}_i{rec['i']}_j{rec['j']}_{ph}"
    ph = "pre" if rec["phase"] == 3 else "post"
    return f"sed_r{rank}_t{rec['TIME_AMPS']}_i{rec['i']}_j{rec['j']}_s{rec['isn']}_{ph}"


def aggregate(dump_dir: str | Path) -> dict:
    """Parse every amps_dump_r*_t*.bin file under dump_dir and merge all
    records into one flat dict of rank-qualified keys -> numpy arrays /
    scalars. Raises ValueError if two source files ever produce the same
    key (which would otherwise silently drop data — see module docstring)."""
    dump_dir = Path(dump_dir)
    files = sorted(dump_dir.glob("amps_dump_r*_t*.bin"))
    if not files:
        raise SystemExit(f"no amps_dump_r*_t*.bin files in {dump_dir}")
    out: dict = {}
    key_src: dict[str, Path] = {}
    for f in files:
        rank, _thread = _parse_fname(f)
        for rec in read_dump_file(f):
            base = _key(rec, rank)
            for name, val in rec.items():
                if name == "kind":
                    continue
                key = f"{base}_{name}"
                if key in out:
                    raise ValueError(
                        f"duplicate key {key!r}: written by both {key_src[key]} and {f}"
                    )
                out[key] = val
                key_src[key] = f
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump_dir", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("amps_ref.npz"))
    args = ap.parse_args()
    out = aggregate(args.dump_dir)
    n_files = len(sorted(args.dump_dir.glob("amps_dump_r*_t*.bin")))
    np.savez_compressed(args.output, **out)
    print(f"wrote {args.output} ({len(out)} arrays from {n_files} files)")


if __name__ == "__main__":
    main()
