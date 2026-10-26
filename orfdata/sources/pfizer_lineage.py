"""Reconcile deposited dog run/channel identities with published tissue pools.

These functions only inspect caller-supplied source facts. The vendor's generic
``Sample`` label is preserved separately from S2's tissue/bridge/empty roles.
No individual animal identity, positive donor or technical repeat is inferred.
"""

import re
from collections import Counter

from ..measurements import Assay, Sample
from .pfizer import SAMPLE_COLUMNS, protein_sample_id

DOG_FILESETS = {f"IQP_1269_Dog_P{i}": f"Plex{i}" for i in range(1, 5)}
CHANNELS = ("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C")
ACCESSION = "MSV000099744"


def native_sample_id(channel):
    return ACCESSION + ":pd:sample:" + channel["sample_id"]


def native_assay_id(channel, run=None):
    if run is None:
        return native_sample_id(channel) + ":aggregate"
    return f"{ACCESSION}:pd:run:{run['workflow_id']}:{run['file_id']}:channel:{channel['channel_id']}"


def lineage_records(lineage, reference):
    """Expose reconciled source pools and run/fraction channels in the common API.

    Aggregate assays carry 44-sample peptide-group arrays. Run/channel assays
    carry PSM/reporter observations. They share biological samples across all
    fractions; technical fractions do not manufacture biological replicates.
    Empty and bridge channels have no tissue or individual-dog membership.
    """
    if reference.taxon_id != 9615:
        raise ValueError("Dog lineage requires a canine reference")
    samples, assays = [], []
    provenance = dict(study=ACCESSION, source_sha256=lineage["source_sha256"])
    for channel in lineage["channels"]:
        identifier = native_sample_id(channel)
        pool = channel["tissue_pool_id"]
        samples.append(Sample(
            identifier, reference, biological_unit_id=ACCESSION + ":pool:" + pool if pool else None,
            tissue=channel["tissue"], condition=channel["role"], pool_members=(),
            attributes=dict(role=channel["role"], pool_membership=channel["pool_membership"],
                            source_channel=channel, donor_id=None, sex=None, biological_strain=None,
                            individual_positive_dogs="not identifiable from a tissue pool"),
            provenance=provenance))
        common = dict(sample_id=identifier, kind="digest_ms", reference=reference,
                      channel=channel["channel"], plex=channel["plex"], bridge="131C",
                      search_reference=reference.key,
                      normalization=dict(method="native PD column; see each observation's unit and provenance",
                                         transformed_by_importer=False, absolute_molecule_counts=False),
                      calibration=dict(method="none to copies/cell"), provenance=provenance)
        summary = channel["summary_sample"]
        assays.append(Assay(
            assay_id=native_assay_id(channel), **common,
            attributes=dict(role=channel["role"], aggregation="source peptide-group channel across fractions",
                            file_set_number=channel["file_set_number"], source_sample_id=channel["sample_id"],
                            protein_summary_sample_id=channel["summary_sample_id"],
                            summary_technical_replicate=summary["TechRepMixture"] if summary else None,
                            instrument_reconciliation=lineage["instrument_reconciliation"])))
        for run in lineage["runs"]:
            if run["file_set_id"] != channel["file_set_id"]:
                continue
            assays.append(Assay(
                assay_id=native_assay_id(channel, run), **common,
                run_id=f"{ACCESSION}:pd:run:{run['workflow_id']}:{run['file_id']}", fraction=str(run["fraction"]),
                attributes=dict(role=channel["role"], workflow_id=run["workflow_id"], file_id=run["file_id"],
                                study_file_id=run["study_file_id"], file_name=run["file_name"],
                                native_instrument=run["instrument"], software_revision=run["software_revision"],
                                technical_replicate_status=run["technical_replicate_status"],
                                instrument_reconciliation=lineage["instrument_reconciliation"])))
    return tuple(samples), tuple(assays)


def _unique(rows, key):
    result = {}
    for row in rows:
        identifier = key(row)
        if identifier in result:
            raise ValueError("Duplicate source lineage identity: " + repr(identifier))
        result[identifier] = row
    return result


def _basename(path):
    return path.replace("\\", "/").rsplit("/", 1)[-1]


def reconcile_dog_lineage(study, panel, summary_samples, workflow_files):
    """Return checked native identities and all 44 published channel roles.

    The native P1..P4 to publication Plex1..Plex4 association uses the named dog
    batches, explicitly recorded as a label-based association. Summary mixtures
    are independently matched by their *entire* channel/tissue profile, never
    their names. All 38 source summary samples must match exactly. The database
    run IDs join the study by StudyFileID, checked against basename and size.
    """
    filesets = _unique(study["file_sets"], lambda r: r["Id"])
    dog = _unique((r for r in filesets.values() if r["Name"] in DOG_FILESETS), lambda r: r["Name"])
    if set(dog) != set(DOG_FILESETS):
        raise ValueError("Expected the four named dog file sets")
    all_samples = _unique(study["samples"], lambda r: r["Id"])
    native_channels = _unique(study["channels"], lambda r: (r["QuanMethodId"], r["Id"]))
    publication = _unique(panel["channels"], lambda r: (r["plex"], r["channel"]))
    expected = {(plex, channel) for plex in DOG_FILESETS.values() for channel in CHANNELS}
    if set(publication) != expected:
        raise ValueError("Published dog panel must contain all four eleven-channel plexes")
    if Counter(r["source_label"] for r in publication.values())["Bridge"] != 4 or any(
            publication[plex, "131C"]["source_label"] != "Bridge" for plex in DOG_FILESETS.values()):
        raise ValueError("Published bridge roles disagree with the inspected source")
    empty = {key for key, row in publication.items() if row["source_label"] == "Empty"}
    if empty != {("Plex1", "131N"), ("Plex4", "131N")}:
        raise ValueError("Published empty channels disagree with the inspected source")

    summaries = _unique(summary_samples, protein_sample_id)
    if len(summaries) != 38:
        raise ValueError("Expected 38 distinct source protein-summary samples")
    mixture_profiles, summary_by_channel = {}, {}
    for identifier, row in summaries.items():
        if set(row) != set(SAMPLE_COLUMNS) or row["sample.type"] != "sample":
            raise ValueError("Unexpected protein-summary sample fields or role")
        channel = row["Channel"].replace("_", "")
        key = row["Mixture"], channel
        if key in summary_by_channel or channel not in CHANNELS:
            raise ValueError("Duplicate or unknown summary mixture/channel")
        summary_by_channel[key] = identifier, row
        mixture_profiles.setdefault(row["Mixture"], set()).add((channel, row["Condition"]))
    if len(mixture_profiles) != 4:
        raise ValueError("Expected four protein-summary mixtures")
    matched_mixtures = {}
    for plex in DOG_FILESETS.values():
        profile = {(channel, row["source_label"]) for (p, channel), row in publication.items()
                   if p == plex and row["source_label"] not in ("Empty", "Bridge")}
        matches = [name for name, value in mixture_profiles.items() if value == profile]
        if len(matches) != 1:
            raise ValueError("No unique complete channel/tissue profile matches " + plex)
        matched_mixtures[plex] = matches[0]
    if len(set(matched_mixtures.values())) != 4:
        raise ValueError("Protein-summary mixtures do not map one-to-one to publication plexes")

    channels, runs, used_samples, expected_files = [], [], set(), {}
    for name, file_set in sorted(dog.items()):
        plex = DOG_FILESETS[name]
        file_set_number = int(file_set["SampleFileSetNumber"])
        native = [s for s in all_samples.values() if s["FileSetId"] == file_set["Id"]]
        by_channel = {}
        for sample in native:
            channel_info = sample["quan_channel"]
            key = channel_info["QuanMethodId"], channel_info["QuanChannelId"]
            if channel_info["QuanMethodId"] != file_set["QuanMethodId"] or key not in native_channels:
                raise ValueError("Study sample quantification method/channel mismatch")
            definition = native_channels[key]
            channel = definition["Name"]
            if channel in by_channel:
                raise ValueError("Duplicate native file-set channel")
            by_channel[channel] = sample
            source = publication[plex, channel]
            label = source["source_label"]
            role = {"Bridge": "bridge", "Empty": "empty"}.get(label, "tissue_pool")
            summary_id, summary = summary_by_channel.get((matched_mixtures[plex], channel), (None, None))
            if (role == "tissue_pool") != (summary is not None):
                raise ValueError("Tissue/control role conflicts with protein-summary coverage")
            if summary_id is not None:
                used_samples.add(summary_id)
            channels.append(dict(
                sample_id=sample["Id"], source_sample_name=sample["Name"],
                source_sample_number=int(sample["SampleNumber"]), vendor_sample_type=sample["SampleType"],
                file_set_id=file_set["Id"], file_set_name=name, file_set_number=file_set_number,
                plex=plex, channel=channel, channel_id=int(definition["Id"]),
                channel_position=int(definition["Position"]), quantification_method_id=int(definition["QuanMethodId"]),
                role=role, source_label=label, tissue=label if role == "tissue_pool" else None,
                tissue_pool_id=summary["BioReplicate"] if summary else None,
                summary_sample_id=summary_id, summary_sample=summary,
                pool_members=[], pool_membership="unresolved" if role != "empty" else "not_applicable",
                publication_cell=source["cell"], publication_sheet=panel["sheet"],
                native_to_publication_basis="matching numbered dog batch labels (P1..P4 / Plex1..Plex4)",
                summary_to_publication_basis="unique complete channel/tissue profile equality"))
        if set(by_channel) != set(CHANNELS):
            raise ValueError("Study dog file set must have eleven distinct channels")
        fractions = _unique(file_set["files"], lambda r: int(r["FractionNumber"]))
        if set(fractions) != set(range(1, 25)):
            raise ValueError("Expected 24 distinct source fractions per dog file set")
        for fraction, file in sorted(fractions.items()):
            study_id = f"F{file_set_number}.{fraction}"
            if study_id in expected_files:
                raise ValueError("Duplicate study file identity")
            expected_files[study_id] = (file_set, file, plex, fraction)
    if used_samples != set(summaries):
        raise ValueError("Not every source protein-summary sample was reconciled")

    raw_files = _unique((r for r in workflow_files if r["FileType"] == "XcaliburRawfile"),
                        lambda r: r["StudyFileID"])
    if set(raw_files) != set(expected_files):
        raise ValueError("Native run coverage differs from the 96 study files")
    for study_id, row in sorted(raw_files.items()):
        file_set, file, plex, fraction = expected_files[study_id]
        size = re.fullmatch(r"(\d+) \[Byte\]", file["FileSize"])
        if (_basename(row["FileName"]) != _basename(file["FileName"]) or size is None
                or int(size[1]) != row["FileSize"]):
            raise ValueError("Native run identity/size disagrees with study metadata: " + study_id)
        runs.append(dict(workflow_id=row["WorkflowID"], file_id=row["FileID"], study_file_id=study_id,
                         file_set_id=file_set["Id"], file_set_name=file_set["Name"],
                         file_set_number=int(file_set["SampleFileSetNumber"]), plex=plex,
                         fraction=fraction, file_name=_basename(row["FileName"]), file_size=row["FileSize"],
                         instrument=row["InstrumentDescription"], software_revision=row["SoftwareRevision"],
                         technical_replicate=None, technical_replicate_status="not independently identified",
                         vendor=row, study_file=file))
    _unique(runs, lambda r: (r["workflow_id"], r["file_id"]))
    return dict(channels=channels, runs=runs,
                counts=dict(channels=len(channels), runs=len(runs),
                            **Counter(c["role"] for c in channels)),
                source_sha256=dict(study=study["source_sha256"], supplementary_workbook=panel["source_sha256"]),
                instrument_reconciliation=dict(
                    status="conflicting source reports; native run metadata supports deposit report",
                    publication=dict(instrument="Orbitrap Fusion Eclipse", citation="https://doi.org/10.1038/s41597-026-07512-9"),
                    deposit=dict(instrument="Orbitrap Fusion Lumos", cv_accession="MS:1002732",
                                 citation="https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD070290"),
                    native_run_counts=dict(Counter(r["instrument"] for r in runs))))
