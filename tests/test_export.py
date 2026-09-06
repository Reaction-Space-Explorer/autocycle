import pathlib
import sys

import pytest

import autocycle.export as export
from autocycle.export import _load_cairosvg, write
from autocycle.io_spec import load_yaml
from autocycle.render import render

SVG = "<svg xmlns='http://www.w3.org/2000/svg' width='10' height='10' viewBox='0 0 10 10'/>"


def test_svg_is_written_verbatim(tmp_path):
    out = write(SVG, tmp_path / "a.svg")
    assert out.read_text() == SVG


def test_unknown_suffix_is_treated_as_svg(tmp_path):
    assert write(SVG, tmp_path / "a.txt").read_text() == SVG


def test_png_and_pdf_when_cairosvg_works(tmp_path):
    try:
        _load_cairosvg()
    except Exception as exc:  # noqa: BLE001 - cairocffi raises OSError without libcairo
        pytest.skip(f"cairosvg unavailable: {exc}")
    svg = render(load_yaml("examples/formose_gain.yaml"))
    png = write(svg, tmp_path / "a.png", width=200)
    pdf = write(svg, tmp_path / "a.pdf")
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert pdf.read_bytes()[:5] == b"%PDF-"


def test_homebrew_libcairo_is_found_without_a_preset_env(monkeypatch):
    if sys.platform != "darwin" or not any(
        (pathlib.Path(d) / "libcairo.2.dylib").exists() for d in export.BREW_LIB
    ):
        pytest.skip("not macOS with a Homebrew libcairo")
    monkeypatch.delenv("DYLD_FALLBACK_LIBRARY_PATH", raising=False)
    assert _load_cairosvg() is not None


def test_the_error_names_the_variable_that_fixes_it(tmp_path, monkeypatch):
    monkeypatch.setattr(export, "_load_cairosvg", lambda: (_ for _ in ()).throw(OSError("no libcairo")))
    monkeypatch.setattr(export.sys, "platform", "darwin")
    with pytest.raises(RuntimeError, match="DYLD_FALLBACK_LIBRARY_PATH"):
        write(SVG, tmp_path / "a.png")
