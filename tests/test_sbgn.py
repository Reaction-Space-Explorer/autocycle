import xml.etree.ElementTree as ET
from collections import Counter

from autocycle.io_spec import load_yaml
from autocycle.sbgn import NS, to_sbgn

SPEC = "examples/canonical/formose_core.yaml"
NSMAP = {"s": NS}


def _root(path=SPEC):
    return ET.fromstring(to_sbgn(load_yaml(path)))


def test_is_well_formed_sbgn_ml():
    r = _root()
    assert r.tag == f"{{{NS}}}sbgn"
    assert r.find("s:map", NSMAP).get("language") == "process description"


def test_molecules_are_simple_chemicals_and_reactions_are_processes():
    cls = Counter(g.get("class") for g in _root().findall(".//s:glyph", NSMAP))
    c = load_yaml(SPEC)
    n_side = sum(len(s.consumes) + len(s.produces) for s in c.steps)
    assert cls["process"] == len(c.steps)
    assert cls["simple chemical"] == len(c.nodes) + n_side


def test_every_ring_step_has_a_consumption_and_a_production_arc():
    arcs = _root().findall(".//s:arc", NSMAP)
    cls = Counter(a.get("class") for a in arcs)
    assert cls["consumption"] >= len(load_yaml(SPEC).steps)
    assert cls["production"] >= len(load_yaml(SPEC).steps)
    assert set(cls) <= {"consumption", "production"}


def test_arcs_reference_declared_glyphs():
    r = _root()
    ids = {g.get("id") for g in r.findall(".//s:glyph", NSMAP)}
    for a in r.findall(".//s:arc", NSMAP):
        assert a.get("source") in ids
        assert a.get("target") in ids


def test_glyphs_carry_a_bbox_and_a_label():
    for g in _root().findall(".//s:glyph", NSMAP):
        assert g.find("s:bbox", NSMAP) is not None
        if g.get("class") == "simple chemical":
            assert g.find("s:label", NSMAP).get("text")


def test_layout_coordinates_are_distinct():
    xs = {
        (g.find("s:bbox", NSMAP).get("x"), g.find("s:bbox", NSMAP).get("y"))
        for g in _root().findall(".//s:glyph", NSMAP)
    }
    assert len(xs) > 4


def test_the_lossy_parts_are_recorded_as_notes():
    """SBGN-PD cannot express a shunt or a gain, so they must not vanish silently."""
    body = _root().find(".//s:notes/s:body", NSMAP)
    assert "no glyph for a shunt" in "".join(body.itertext())
    assert "autocatalyst at ring position 0" in "".join(body.itertext())
