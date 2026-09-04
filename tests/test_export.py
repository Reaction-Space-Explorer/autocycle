import pytest

from autocycle.export import write
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
        import cairosvg  # noqa: F401
    except Exception as exc:  # noqa: BLE001 - cairocffi raises OSError without libcairo
        pytest.skip(f"cairosvg unavailable: {exc}")
    svg = render(load_yaml("examples/formose_gain.yaml"))
    png = write(svg, tmp_path / "a.png", width=200)
    pdf = write(svg, tmp_path / "a.pdf")
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert pdf.read_bytes()[:5] == b"%PDF-"
