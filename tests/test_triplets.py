import pytest

from autocycle.ingest import find_cycles, to_cycle
from autocycle.io_triplets import gain_step, orient, read_triplets
from autocycle.spec import SpecError
from autocycle.verify import AUTOCATALYTIC, verify

# the formose core as a stoichiometric matrix: glycolaldehyde, glyceraldehyde,
# dihydroxyacetone, erythrose, with two formaldehyde in and two glycolaldehyde out
CORE = """\
r1\tOCC=O\t-1\tAldol addition
r1\tC=O\t-1\tAldol addition
r1\tOCC(O)C=O\t1\tAldol addition
r2\tOCC(O)C=O\t-1\tKeto-enol
r2\tOCC(=O)CO\t1\tKeto-enol
r3\tOCC(=O)CO\t-1\tAldol addition
r3\tC=O\t-1\tAldol addition
r3\tOCC(O)C(O)C=O\t1\tAldol addition
r4\tOCC(O)C(O)C=O\t-1\tRetro-aldol
r4\tOCC=O\t2\tRetro-aldol
"""


@pytest.fixture
def core(tmp_path):
    p = tmp_path / "core.tsv"
    p.write_text(CORE)
    return read_triplets(p)


def test_a_cycle_is_oriented_on_the_species_it_amplifies(core):
    """A cycle search returns an arbitrary rotation, so the seed is read from the
    coefficients; otherwise the verdict depends on where the search happened to start."""
    ring = next(c for c in find_cycles(core, 3, 5) if len(c) == 4)
    turned, gain = orient(core, ring)
    assert turned[0] == "O=CCO"
    assert gain == len(ring) - 1
    assert gain_step(core, turned) == gain


def test_the_matrix_reads_back_as_an_autocatalytic_cycle(core):
    ring, gain = orient(core, next(c for c in find_cycles(core, 3, 5) if len(c) == 4))
    cyc = to_cycle(core, ring, gain_at=gain)
    v = verify(cyc)
    assert v.status == AUTOCATALYTIC
    assert v.seed_yield == 2.0


def test_multiplicity_may_be_a_coefficient_or_repeated_rows(tmp_path):
    a = tmp_path / "a.tsv"
    a.write_text(CORE)
    b = tmp_path / "b.tsv"
    b.write_text(CORE.replace("r4\tOCC=O\t2\tRetro-aldol\n",
                              "r4\tOCC=O\t1\tRetro-aldol\nr4\tOCC=O\t1\tRetro-aldol\n"))
    ga, gb = read_triplets(a), read_triplets(b)
    assert sorted(find_cycles(ga, 3, 5)) == sorted(find_cycles(gb, 3, 5))
    ring, gain = orient(gb, next(c for c in find_cycles(gb, 3, 5) if len(c) == 4))
    assert verify(to_cycle(gb, ring, gain_at=gain)).seed_yield == 2.0


def test_a_species_consumed_and_produced_equally_makes_no_edge(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("r1\tOCC=O\t-1\tx\nr1\tOCC=O\t1\tx\nr1\tC=O\t-1\tx\nr1\tOCC(O)C=O\t1\tx\n")
    g = read_triplets(p)
    assert not g.has_edge("O=CCO", "O=CC(O)CO")


def test_a_header_is_used_when_present(tmp_path):
    p = tmp_path / "h.csv"
    p.write_text("reaction,species,coefficient,rule\n" + CORE.replace("\t", ","))
    assert read_triplets(p).number_of_nodes() == 5


def test_a_non_numeric_coefficient_is_named(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("r1\tOCC=O\tmuch\tx\n")
    with pytest.raises(SpecError, match="non-numeric coefficient"):
        read_triplets(p)


def test_an_empty_table_is_named(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text("")
    with pytest.raises(SpecError, match="no rows"):
        read_triplets(p)
