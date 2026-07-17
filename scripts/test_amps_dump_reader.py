#!/usr/bin/env python3
"""Round-trip test for amps_dump_reader: build synthetic records byte-identically
to the Fortran writers (little- and big-endian variants), parse and compare,
and exercise the multi-rank aggregation path (rank-qualified keys, duplicate
guard). Run: python scripts/test_amps_dump_reader.py"""
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from amps_dump_reader import MAGIC_MICRO, MAGIC_SED, aggregate, read_dump_file


def w_i0(b, v, bo="<"):
    return b + struct.pack(f"{bo}i", v)


def w_i1(b, a, bo="<"):
    a = np.asarray(a, dtype=f"{bo}i4")
    return b + struct.pack(f"{bo}i", a.size) + a.tobytes()


def w_r0(b, v, bo="<"):
    return b + struct.pack(f"{bo}d", v)


def w_rn(b, a, bo="<"):
    """Array with per-rank int32 size prefixes, Fortran order."""
    a = np.asarray(a, dtype=f"{bo}f8")
    for s in a.shape:
        b += struct.pack(f"{bo}i", s)
    return b + a.tobytes(order="F")


def build_micro(phase, nmic=3, npr=6, nbr=4, ncr=1, npi=18, nbi=2, nci=1,
                npa=5, nba=1, nca=4, mxnbin=4, bo="<"):
    rng = np.random.default_rng(phase)
    b = b""
    for v in (MAGIC_MICRO, 1, phase, 42, 3, 4, 1, nmic,
              npr, nbr, ncr, npi, nbi, nci, npa, nba, nca, mxnbin,
              1, 12345, 0, 7, 99):
        b = w_i0(b, v, bo)
    b = w_r0(b, 1.0, bo)
    kmicvm = np.arange(2, 2 + nmic)
    b = w_i1(b, kmicvm, bo)
    # dt/kmicvm are tracked here (not just consumed silently) so the
    # round-trip loop below asserts them directly, same as every other field.
    fields = {"dt": 1.0, "kmicvm": kmicvm}
    for name in ("qcvm", "v3v", "qvvm", "moist_denvm", "ptotvm", "tvm", "wbvm",
                 "trpv_thil", "trpv_qtp"):
        fields[name] = rng.uniform(size=nmic)
        b = w_rn(b, fields[name], bo)
    fields["qrpvm"] = rng.uniform(size=(npr, nbr, ncr, nmic))
    fields["qipvm"] = rng.uniform(size=(npi, nbi, nci, nmic))
    fields["qapvm"] = rng.uniform(size=(npa, nba, nca, nmic))
    for name in ("qrpvm", "qipvm", "qapvm"):
        b = w_rn(b, fields[name], bo)
    if phase == 2:
        fields["dmtendlm"] = rng.uniform(size=(10, 2, nmic))
        fields["dcontendlm"] = rng.uniform(size=(10, 2, nmic))
        fields["dbintendlm"] = rng.uniform(size=(7, 2, mxnbin, nmic))
        for name in ("dmtendlm", "dcontendlm", "dbintendlm"):
            b = w_rn(b, fields[name], bo)
    return b, fields


def build_sed(phase, isn=0, np_=6, nb=4, nc=1, nzh=8, bo="<"):
    rng = np.random.default_rng(100 + phase)
    b = b""
    for v in (MAGIC_SED, 1, phase, 42, 3, 4, 1, isn, 1, np_, nb, nc, 2, 5, 2, 6):
        b = w_i0(b, v, bo)
    b = w_r0(b, 1.0, bo)
    # Non-constant so a wrong reshape order ("C" instead of "F") would be
    # caught when nc > 1, not just a wrong shape.
    k1b_flat = np.arange(10, 10 + nb * nc)
    k2b_flat = np.arange(50, 50 + nb * nc)
    b = w_i1(b, k1b_flat, bo)
    b = w_i1(b, k2b_flat, bo)
    fields = {
        "dt": 1.0,
        "k1b": k1b_flat.reshape((nb, nc), order="F"),
        "k2b": k2b_flat.reshape((nb, nc), order="F"),
    }
    fields["qpv"] = rng.uniform(size=(np_, nb, nc, nzh))
    b = w_rn(b, fields["qpv"], bo)
    r1names = ("q_this", "q_other", "qcv", "qtp", "moist_denv", "thetav", "qvv",
               "tv", "dens_col", "momz_col", "u_col", "v_col", "cz_col", "fz_col",
               "dzzmv", "dzvmv")
    for name in r1names:
        fields[name] = rng.uniform(size=nzh)
        b = w_rn(b, fields[name], bo)
    fields["mmass"] = rng.uniform(size=(nb, nc, nzh))
    b = w_rn(b, fields["mmass"], bo)
    for name in ("den_t", "momz_t", "rhou_t", "rhov_t", "rhoe_t"):
        fields[name] = rng.uniform(size=nzh)
        b = w_rn(b, fields[name], bo)
    b = w_r0(b, 3.5, bo)
    fields["sflx"] = 3.5
    return b, fields


def test_round_trip():
    """Original single-rank, mixed-record round trip, plus direct shape
    checks on the fields that used to be read but never independently
    asserted (dt, kmicvm, k1b/k2b)."""
    blob1, f1 = build_micro(1)
    blob2, f2 = build_micro(2)
    blob3, f3 = build_sed(3, nc=2)  # nc>1 exercises the (nb, nc) F-order reshape
    blob4, f4 = build_sed(4, isn=1, nc=2)
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

    # Direct asserts (previously only indirectly exercised, if at all).
    assert recs[0]["dt"] == 1.0
    assert np.array_equal(recs[0]["kmicvm"], np.arange(2, 2 + 3))
    assert recs[1]["k1b"].shape == (4, 2)
    assert recs[1]["k2b"].shape == (4, 2)
    assert np.array_equal(recs[1]["k1b"], f3["k1b"])
    assert np.array_equal(recs[1]["k2b"], f3["k2b"])
    print("OK: 4 records round-tripped (dt/kmicvm/k1b/k2b asserted directly)")


def test_byteswap():
    """A big-endian file (as produced by SCALE cluster builds compiled with
    -fconvert=big-endian / -convert big_endian) must parse identically to
    the little-endian one, via auto-detection off the record magic."""
    blob_le, fields_le = build_micro(1, bo="<")
    blob_be, fields_be = build_micro(1, bo=">")
    assert all(np.allclose(fields_le[k], fields_be[k]) for k in fields_le), (
        "test bug: le/be field builders diverged"
    )
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "amps_dump_r000002_t001.bin"
        p.write_bytes(blob_be)
        recs = read_dump_file(p)
    assert len(recs) == 1, f"expected 1 record, got {len(recs)}"
    rec = recs[0]
    for name, val in fields_be.items():
        assert np.allclose(rec[name], val), f"byteswapped field {name} mismatch"
    print("OK: big-endian record auto-detected and parsed identically to little-endian")


def _n_data_keys(rec: dict) -> int:
    return sum(1 for k in rec if k != "kind")


def test_aggregate_two_ranks():
    """Two ranks dumping the SAME local (i, j) box (the normal case — indices
    are local per rank) must both survive aggregation under distinct,
    rank-qualified keys, not collide into one."""
    blob, fields = build_micro(1)
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "amps_dump_r000000_t001.bin").write_bytes(blob)
        (d / "amps_dump_r000001_t001.bin").write_bytes(blob)
        out = aggregate(d)
        n_per_record = _n_data_keys(read_dump_file(d / "amps_dump_r000000_t001.bin")[0])

    assert len(out) == 2 * n_per_record, (
        f"expected {2 * n_per_record} keys from 2 ranks, got {len(out)}"
    )
    key_r0 = "micro_r0_t42_i3_j4_pre_dt"
    key_r1 = "micro_r1_t42_i3_j4_pre_dt"
    assert key_r0 in out, f"{key_r0!r} missing from aggregate output"
    assert key_r1 in out, f"{key_r1!r} missing from aggregate output"
    assert out[key_r0] == fields["dt"] == out[key_r1]
    print(f"OK: aggregate() kept both ranks' data ({len(out)} keys, no collision)")


def test_duplicate_key_guard():
    """Two files that resolve to the same rank and produce the same record
    key (here: two different thread-files under rank 0 carrying the same
    (TIME_AMPS, i, j) column — not supposed to happen in a real run since a
    column belongs to exactly one thread, but it is exactly the condition
    the guard exists to catch) must raise, not silently drop data."""
    blob, _ = build_micro(1)
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "amps_dump_r000000_t001.bin").write_bytes(blob)
        (d / "amps_dump_r000000_t002.bin").write_bytes(blob)
        try:
            aggregate(d)
        except ValueError as e:
            msg = str(e)
            assert "duplicate key" in msg, f"unexpected ValueError message: {msg}"
            assert "t001" in msg and "t002" in msg, (
                f"guard message should name both source files: {msg}"
            )
        else:
            raise AssertionError("expected ValueError for duplicate key, none was raised")
    print("OK: duplicate-key guard raised ValueError naming both source files")


def main():
    test_round_trip()
    test_byteswap()
    test_aggregate_two_ranks()
    test_duplicate_key_guard()
    print("OK: all amps_dump_reader checks passed "
          "(round-trip, big-endian auto-detect, 2-rank aggregate, duplicate-key guard)")


if __name__ == "__main__":
    main()
