import pytest

from autocycle.encode import GREY, dg_colour, dg_span, widths


def test_widths_monotonic_and_bounded():
    w = widths([1.0, 2.0, 4.0], lo=0.1, hi=0.5)
    assert w == sorted(w)
    assert w[0] == pytest.approx(0.1)
    assert w[-1] == pytest.approx(0.5)


def test_zero_magnitude_stays_zero():
    assert widths([0.0, 1.0, 2.0])[0] == 0.0


def test_all_equal_maps_to_midpoint():
    assert widths([3.0, 3.0], lo=0.1, hi=0.5) == [pytest.approx(0.3)] * 2


def test_log_mode_compresses_decades():
    lin = widths([1.0, 10.0, 100.0], mode="linear")
    log = widths([1.0, 10.0, 100.0], mode="log")
    mid = (lin[0] + lin[2]) / 2
    # linear squashes the middle decade toward the bottom; log centres it
    assert lin[1] < 0.6 * mid
    assert log[1] == pytest.approx((log[0] + log[2]) / 2)


def test_multiples_mode_preserves_ratios():
    w = widths([1.0, 2.0], mode="multiples", hi=0.4)
    assert w[1] / w[0] == pytest.approx(2.0)


def test_unknown_mode_rejected():
    with pytest.raises(ValueError, match="unknown mode"):
        widths([1.0], mode="sqrt")


def test_dg_span_ignores_unknowns():
    assert dg_span([-5.0, None, 12.0]) == 12.0
    assert dg_span([None, None]) == 1.0


def test_dg_colour_unknown_is_grey():
    assert dg_colour(None, 10.0) == GREY


def test_dg_colour_diverges_around_zero():
    blue, mid, red = dg_colour(-10, 10), dg_colour(0, 10), dg_colour(10, 10)
    assert blue != red
    # blue end is bluest, red end is reddest
    assert int(blue[5:7], 16) > int(blue[1:3], 16)
    assert int(red[1:3], 16) > int(red[5:7], 16)
    assert mid != blue and mid != red
