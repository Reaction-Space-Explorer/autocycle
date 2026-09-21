"""Choices the scripts share, stated once.

The food set is a modelling decision rather than a parameter: a count taken under
a different one is a different quantity, so it is written here and nowhere else.
Each network is read at the deepest generation the deposition carries, which is
why the generation is part of the name. Scripts that read a different set, or a
different sweep of generations, say so themselves.
"""

FOOD = {"O", "C=O", "C(=O)=O", "N"}
BASE = {"O", "C=O"}                       # water and formaldehyde, before CO2 and NH3

#                       reaction table under RELS        energy-table stem
NETS = {
    "glucose G5":         ("Glucose/GlucoseRels_5.tsv", "Glucose_G5"),
    "glucose+ammonia G4": ("GlucoseAmm/GlucoseAmmRels_4.tsv", "GlucoseAmm_G4"),
    "formose G6":         ("Formose/FormoseRels_6.tsv", "Formose_G6"),
    "formose+ammonia G4": ("FormoseAmm/FormoseAmmRels_4.tsv", "FormoseAmm_G4"),
    "pyruvic acid G6":    ("PyruvicAcid/PyruvicAcidRels_6.tsv", "PyruvicAcid_G6"),
}
