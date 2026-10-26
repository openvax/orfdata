"""Read the pinned dog study/S2/summary sources and export checked pool lineage."""

import argparse
import hashlib
import json
import posixpath
import sqlite3
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from orfdata.sources.pfizer import SAMPLE_COLUMNS, _rows, protein_sample_id
from orfdata.sources.pfizer_lineage import reconcile_dog_lineage

STUDY_SHA256 = "fd0a494cd0a73d120045e092668f7ce02cc23ec57722cc48d309e57ca6394bb8"
PANEL_SHA256 = "a95b7a2514679ba6f7283573059f27e6182b5620f6d0beead601e47daf85e740"
SUMMARY_SHA256 = "ae2953c3fe3c31bc5fa06a9a0b222815376faea25489e6626af3736b7f4cbef2"
SCHEMA_SHA256 = "4748720e8b75b1399c7dee6b40bd9a24ddff8a734e39c07a01835b1fab5ed9c1"
SHEET = "TableS2_tissues_proteomics"
SML = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"
DOCREL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def file_sha256(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def study_metadata(path):
    path = Path(path)
    if path.stat().st_size > 8 * 1024 * 1024:
        raise ValueError("Study metadata exceeds the explicit 8 MiB limit")
    data = path.read_bytes()
    if b"<!DOCTYPE" in data or b"<!ENTITY" in data:
        raise ValueError("Study XML must not declare entities")
    definition = ET.fromstring(data).find("StudyDefinition")
    if definition is None:
        raise ValueError("Missing native StudyDefinition")
    return dict(
        source_sha256=hashlib.sha256(data).hexdigest(),
        file_sets=[dict(e.attrib, files=[dict(f.attrib) for f in e.findall("Files/File")])
                   for e in definition.findall("FileSets/FileSet")],
        channels=[dict(e.attrib) for e in definition.findall("QuanMethods/QuanMethod/QuanChannels/QuanChannel")],
        samples=[dict(e.attrib, quan_channel=dict(e.find("QuanChannelInformation").attrib))
                 for e in definition.findall("Samples/Sample")])


def dog_panel(path):
    """Read literal S2 cells M1:Q13 without changing or recalculating the workbook."""
    path = Path(path)
    if path.stat().st_size > 2 * 1024 * 1024:
        raise ValueError("Supplement exceeds the explicit 2 MiB compressed limit")
    with zipfile.ZipFile(path) as archive:
        if sum(info.file_size for info in archive.infolist()) > 16 * 1024 * 1024:
            raise ValueError("Supplement exceeds the explicit 16 MiB expanded limit")
        if len(archive.namelist()) != len(set(archive.namelist())):
            raise ValueError("Duplicate workbook members")

        def xml(name):
            data = archive.read(name)
            if b"<!DOCTYPE" in data or b"<!ENTITY" in data:
                raise ValueError("Workbook XML must not declare entities")
            return ET.fromstring(data)

        sheets = [s for s in xml("xl/workbook.xml").findall(f"{SML}sheets/{SML}sheet") if s.get("name") == SHEET]
        if len(sheets) != 1:
            raise ValueError("Expected one published S2 worksheet")
        relation_id = sheets[0].get(DOCREL + "id")
        relations = [r for r in xml("xl/_rels/workbook.xml.rels").findall(REL + "Relationship")
                     if r.get("Id") == relation_id and r.get("TargetMode") != "External"]
        if len(relations) != 1:
            raise ValueError("Missing internal S2 worksheet relationship")
        member = posixpath.normpath(posixpath.join("xl", relations[0].get("Target")))
        if not member.startswith("xl/"):
            raise ValueError("Worksheet relationship escapes the workbook")
        strings = ["".join(e.itertext()) for e in xml("xl/sharedStrings.xml").findall(SML + "si")]
        cells = {}
        wanted = {f"{column}{row}" for column in "MNOPQ" for row in range(1, 14)}
        for cell in xml(member).iter(SML + "c"):
            address = cell.get("r")
            if address not in wanted:
                continue
            if address in cells or cell.find(SML + "f") is not None:
                raise ValueError("Expected unique literal published panel cells")
            value = cell.find(SML + "v")
            if value is None:
                continue
            if cell.get("t") == "s":
                cells[address] = strings[int(value.text)]
            elif cell.get("t") in (None, "n"):
                # Channel 126 is a numeric cell in the original workbook.
                cells[address] = str(int(value.text))
            else:
                raise ValueError("Unexpected published panel cell type")
    if cells.get("M1") != "Dog" or [cells.get(c + "1") for c in "NOPQ"] != [f"Plex{i}" for i in range(1, 5)]:
        raise ValueError("Published dog panel layout has changed")
    rows = [dict(plex=cells[column + "1"], channel=cells[f"M{row}"],
                 source_label=cells[f"{column}{row}"], cell=f"{column}{row}")
            for column in "NOPQ" for row in range(2, 13)]
    return dict(source_sha256=file_sha256(path), sheet=SHEET, range="M1:Q13", channels=rows,
                searched_fasta_note=cells["M13"], cells=cells)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("study", "supplement", "protein_summary", "native_export", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    for path, digest in ((args.study, STUDY_SHA256), (args.supplement, PANEL_SHA256),
                         (args.protein_summary, SUMMARY_SHA256)):
        if file_sha256(path) != digest:
            raise ValueError("Pinned source hash mismatch: " + str(path))
    db = sqlite3.connect(args.native_export.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        configuration = json.loads(db.execute("SELECT value FROM _export_metadata WHERE key='configuration'").fetchone()[0])
        if configuration["plan"]["source_schema_sha256"] != SCHEMA_SHA256:
            raise ValueError("Unexpected native extraction schema")
        if db.execute("SELECT complete,rows FROM _export_progress WHERE name='WorkflowInputFiles'").fetchone() != (1, 97):
            raise ValueError("WorkflowInputFiles must be completely extracted before lineage reconciliation")
        db.row_factory = sqlite3.Row
        files = [dict(row) for row in db.execute("SELECT * FROM WorkflowInputFiles ORDER BY _source_rowid")]
    finally:
        db.close()
    samples = {}
    for _, row in _rows(args.protein_summary):
        samples[protein_sample_id(row)] = {key: row[key] for key in SAMPLE_COLUMNS}
    result = reconcile_dog_lineage(study_metadata(args.study), dog_panel(args.supplement), samples.values(), files)
    result["source_sha256"]["protein_summary"] = SUMMARY_SHA256
    result["native_export_configuration"] = configuration
    result["tool_sha256"] = file_sha256(__file__)
    result["normalizer_sha256"] = file_sha256(Path(__file__).parent.parent / "orfdata/sources/pfizer_lineage.py")
    data = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    if len(data) > 2 * 1024 * 1024:
        raise ValueError("Lineage output exceeds its explicit 2 MiB limit")
    with args.output.open("xb") as handle:
        handle.write(data)
    print(json.dumps(dict(counts=result["counts"], output_bytes=len(data), sha256=hashlib.sha256(data).hexdigest())))


if __name__ == "__main__":
    main()
