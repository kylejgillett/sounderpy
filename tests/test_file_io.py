"""Tests for the current SounderPy file I/O API.

This replaces the stale test that referenced ``read_sharppy_file``.  The public
reader in the current package is ``from_sharppy``.
"""

from __future__ import annotations

import numpy as np
import pytest

import sounderpy as spy


REQUIRED_KEYS = {"p", "z", "T", "Td", "u", "v", "wd", "ws", "site_info"}


def test_sharppy_round_trip(clean_profile, tmp_path):
    path = tmp_path / "round_trip.snd"

    spy.to_file("sharppy", clean_profile, filename=str(path))
    assert path.exists()
    assert path.stat().st_size > 0

    loaded = spy.from_sharppy(path)

    assert REQUIRED_KEYS.issubset(loaded)
    assert len(loaded["p"]) == len(clean_profile["p"])
    assert loaded["site_info"]["site-id"] == "TEST"
    assert loaded["site_info"]["valid-time"] == [2024, 5, 21, 0]

    np.testing.assert_allclose(
        loaded["p"].to("hPa").magnitude,
        clean_profile["p"].to("hPa").magnitude,
        rtol=0,
        atol=1e-5,
    )
    np.testing.assert_allclose(
        loaded["T"].to("degC").magnitude,
        clean_profile["T"].to("degC").magnitude,
        rtol=0,
        atol=1e-5,
    )


def test_csv_export_contains_expected_columns(clean_profile, tmp_path):
    path = tmp_path / "profile.csv"
    spy.to_file("csv", clean_profile, filename=str(path))

    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    assert first_line == "p,z,T,Td,u,v"


def test_from_sharppy_rejects_empty_file(tmp_path):
    path = tmp_path / "empty.snd"
    path.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="empty"):
        spy.from_sharppy(path)


def test_from_sharppy_rejects_non_sharppy_file(tmp_path):
    path = tmp_path / "not_sharppy.txt"
    path.write_text("this is not a sounding\n", encoding="utf-8")

    with pytest.raises(ValueError, match="does not appear to be a SHARPpy"):
        spy.from_sharppy(path)
