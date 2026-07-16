#!/usr/bin/env python3
"""Round-trip test for amps_dump_reader: build synthetic records byte-identically
to the Fortran writers, then parse and compare. Run: python scripts/test_amps_dump_reader.py"""
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from amps_dump_reader import MAGIC_MICRO, MAGIC_SED, read_dump_file


def w_i0(b, v):
    b += struct.pack("<i", v)
    return b


def w_i1(b, a):
    a = np.asarray(a, dtype="<i4")
    return b + struct.pack("<i", a.size) + a.tobytes()


def w_r0(b, v):
    return b + struct.pack("<d", v)


def w_rn(b, a):
    """Array with per-rank int32 size prefixes, Fortran order."""
    a = np.asarray(a, dtype="<f8")
    for s in a.shape:
        b += struct.pack("<i", s)
    return b + a.tobytes(order="F")


def build_micro(phase, nmic=3, npr=6, nbr=4, ncr=1, npi=18, nbi=2, nci=1,
                npa=5, nba=1, nca=4, mxnbin=4):
    rng = np.random.default_rng(phase)
    b = b""
    for v in (MAGIC_MICRO, 1, phase, 42, 3, 4, 1, nmic,
              npr, nbr, ncr, npi, nbi, nci, npa, nba, nca, mxnbin,
              1, 12345, 0, 7, 99):
        b = w_i0(b, v)
    b = w_r0(b, 1.0)
    b = w_i1(b, np.arange(2, 2 + nmic))
    fields = {}
    for name in ("qcvm", "v3v", "qvvm", "moist_denvm", "ptotvm", "tvm", "wbvm",
                 "trpv_thil", "trpv_qtp"):
        fields[name] = rng.uniform(size=nmic)
        b = w_rn(b, fields[name])
    fields["qrpvm"] = rng.uniform(size=(npr, nbr, ncr, nmic))
    fields["qipvm"] = rng.uniform(size=(npi, nbi, nci, nmic))
    fields["qapvm"] = rng.uniform(size=(npa, nba, nca, nmic))
    for name in ("qrpvm", "qipvm", "qapvm"):
        b = w_rn(b, fields[name])
    if phase == 2:
        fields["dmtendlm"] = rng.uniform(size=(10, 2, nmic))
        fields["dcontendlm"] = rng.uniform(size=(10, 2, nmic))
        fields["dbintendlm"] = rng.uniform(size=(7, 2, mxnbin, nmic))
        for name in ("dmtendlm", "dcontendlm", "dbintendlm"):
            b = w_rn(b, fields[name])
    return b, fields


def build_sed(phase, isn=0, np_=6, nb=4, nc=1, nzh=8):
    rng = np.random.default_rng(100 + phase)
    b = b""
    for v in (MAGIC_SED, 1, phase, 42, 3, 4, 1, isn, 1, np_, nb, nc, 2, 5, 2, 6):
        b = w_i0(b, v)
    b = w_r0(b, 1.0)
    b = w_i1(b, np.full(nb * nc, 2))
    b = w_i1(b, np.full(nb * nc, 5))
    fields = {"qpv": rng.uniform(size=(np_, nb, nc, nzh))}
    b = w_rn(b, fields["qpv"])
    r1names = ("q_this", "q_other", "qcv", "qtp", "moist_denv", "thetav", "qvv",
               "tv", "dens_col", "momz_col", "u_col", "v_col", "cz_col", "fz_col",
               "dzzmv", "dzvmv")
    for name in r1names:
        fields[name] = rng.uniform(size=nzh)
        b = w_rn(b, fields[name])
    fields["mmass"] = rng.uniform(size=(nb, nc, nzh))
    b = w_rn(b, fields["mmass"])
    for name in ("den_t", "momz_t", "rhou_t", "rhov_t", "rhoe_t"):
        fields[name] = rng.uniform(size=nzh)
        b = w_rn(b, fields[name])
    b = w_r0(b, 3.5)
    fields["sflx"] = 3.5
    return b, fields


def main():
    blob1, f1 = build_micro(1)
    blob2, f2 = build_micro(2)
    blob3, f3 = build_sed(3)
    blob4, f4 = build_sed(4, isn=1)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "amps_dump_r000000_t001.bin"
        p.write_bytes(blob1 + blob3 + blob4 + blob2)
        recs = read_dump_file(p)
    assert len(recs) == 4, f"expected 4 records, got {len(recs)}"
    assert [r["kind"] for r in recs] == ["micro", "sed", "sed", "micro"]
    assert recs[0]["phase"] == 1 and recs[3]["phase"] == 2
    assert recs[0]["TIME_AMPS"] == 42 and recs[0]["i"] == 3 and recs[0]["j"] == 4
    for expected, rec in ((f1, recs[0]), (f2, recs[3]), (f3, recs[1]), (f4, recs[2])):
        for name, val in expected.items():
            got = rec[name]
            assert np.allclose(got, val), f"{rec['kind']} field {name} mismatch"
    assert recs[1]["isn"] == 0 and recs[2]["isn"] == 1
    print("OK: 4 records round-tripped")


if __name__ == "__main__":
    main()
