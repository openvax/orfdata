"""Retention-time standard peptides that labs spike into LC-MS runs.

They are synthetic, so a search database that lacks them can assign their (often very
intense) spectra to whatever similar sequence it has, including candidate ncORFs.
build_db.py adds each one to the canonical class of the search database, and
aggregate.py counts their PSMs per run.
"""

# Biognosys iRT kit (unlabelled); the same list as orfdata.sources.kaabinejadian.IRT_PEPTIDES.
IRT = (
    "LGGNEQVTR", "YILAGVENSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR", "GAGSSEPVTGLDAK",
    "TPVISGGPYEYR", "VEATFGVDESNAK", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
    "LFLQFGAQGSPFLK",
)

# Thermo Pierce Retention Time Calibration (PRTC) mixture. The kit's peptides carry heavy
# C-terminal lysine or arginine, which the Sage parameters here don't search for; see README.
PRTC = (
    "SSAAPPPPPR", "GISNEGQNASIK", "HVLTSIGEK", "DIPVPKPK", "IGDYAGIK", "TASEFDSAIAQDK",
    "SAAGAFGPELSR", "ELGQSGVDTYLQTK", "GLILVGGYGTR", "GILFVGSGVSGGEEGAR", "SFANQPLEVVYSK",
    "LTILEELR", "NGFILDGFPR", "ELASGLSFPVGFK", "LSSEAPALFQFDLK",
)

STANDARDS = tuple(("iRT", p) for p in IRT) + tuple(("PRTC", p) for p in PRTC)
_IL = {p.replace("I", "L"): (kit, p) for kit, p in STANDARDS}


def standard_of(peptide):
    """(kit, standard) if the unmodified peptide is a standard, treating I and L as equal."""
    return _IL.get(peptide.replace("I", "L"))


def contained_standards(sequence):
    """Standards that occur in a protein sequence, treating I and L as equal."""
    il = sequence.replace("I", "L")
    return [p for p_il, (_, p) in _IL.items() if p_il in il]
