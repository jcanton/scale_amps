#!/usr/bin/env python3
"""Parse AMPS reference-data dump files written by the instrumented
scale_atmos_phy_mp_amps.F90 (branch cloudlab_port).

Binary layout v1: sequence of records; ints int32 LE, reals float64 LE,
arrays prefixed by one int32 size per rank, Fortran (column-major) order.
Record types: micro (magic 1095586131) and sed (magic 1095586132) — field
order defined in the M0 plan Tasks 2-3 and mirrored here.

Usage: python amps_dump_reader.py DUMP_DIR -o out.npz
"""
from __future__ import annotations

import argparse
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
    def __init__(self, buf: bytes):
        self.buf = buf
        self.pos = 0

    def eof(self) -> bool:
        return self.pos >= len(self.buf)

    def i4(self, n: int = 1):
        out = np.frombuffer(self.buf, dtype="<i4", count=n, offset=self.pos)
        self.pos += 4 * n
        return int(out[0]) if n == 1 else out.copy()

    def f8(self, n: int = 1):
        out = np.frombuffer(self.buf, dtype="<f8", count=n, offset=self.pos)
        self.pos += 8 * n
        return float(out[0]) if n == 1 else out.copy()

    def i1_arr(self) -> np.ndarray:
        n = self.i4()
        return self.i4(n) if n > 0 else np.empty(0, dtype="<i4")

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
    rec["k1b"] = c.i1_arr()
    rec["k2b"] = c.i1_arr()
    rec["qpv"] = c.rn_arr(4)
    for name in SED_R1:
        rec[name] = c.rn_arr(1)
    rec["mmass"] = c.rn_arr(3)
    for name in SED_R1_TAIL:
        rec[name] = c.rn_arr(1)
    rec["sflx"] = c.f8()
    return rec


def read_dump_file(path: str | Path) -> list[dict]:
    c = _Cursor(Path(path).read_bytes())
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


def _key(rec: dict) -> str:
    if rec["kind"] == "micro":
        ph = "pre" if rec["phase"] == 1 else "post"
        return f"micro_t{rec['TIME_AMPS']}_i{rec['i']}_j{rec['j']}_{ph}"
    ph = "pre" if rec["phase"] == 3 else "post"
    return f"sed_t{rec['TIME_AMPS']}_i{rec['i']}_j{rec['j']}_s{rec['isn']}_{ph}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump_dir", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("amps_ref.npz"))
    args = ap.parse_args()
    out: dict = {}
    files = sorted(args.dump_dir.glob("amps_dump_r*_t*.bin"))
    if not files:
        raise SystemExit(f"no amps_dump_r*_t*.bin files in {args.dump_dir}")
    for f in files:
        for rec in read_dump_file(f):
            base = _key(rec)
            for name, val in rec.items():
                if name == "kind":
                    continue
                out[f"{base}_{name}"] = val
    np.savez_compressed(args.output, **out)
    print(f"wrote {args.output} ({len(out)} arrays from {len(files)} files)")


if __name__ == "__main__":
    main()
