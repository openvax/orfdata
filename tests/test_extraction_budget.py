"""Transfer accounting remains bounded across retries, crashes and restarts."""

import json

import pytest

from tools.extraction_budget import LimitExceeded, ScratchBudget, TransferLedger


def test_reservations_survive_interruption_and_partial_attempts_are_counted(tmp_path):
    path = tmp_path / "transfers.sqlite"
    ledger = TransferLedger(path, {"source": "pinned"})
    interrupted = ledger.reserve(100, 20, byte_limit=40, request_limit=3)
    ledger.close()  # Same durable state as process death after the reservation commit.
    resumed = TransferLedger(path, {"source": "pinned"})
    assert resumed.summary() == dict(charged_payload_bytes=20, verified_received_payload_bytes=0,
                                     range_requests=1, unsettled_requests=1, uncertain_reserved_bytes=20)
    partial = resumed.reserve(100, 20, byte_limit=40, request_limit=3)
    resumed.settle(partial, 5, status="failed")
    with pytest.raises(LimitExceeded, match="payload"):
        resumed.reserve(100, 20, byte_limit=40, request_limit=3)
    retry = resumed.reserve(100, 15, byte_limit=40, request_limit=3)
    resumed.settle(retry, 15, status="complete")
    with pytest.raises(LimitExceeded, match="request"):
        resumed.reserve(100, 1, byte_limit=100, request_limit=3)
    assert resumed.summary()["charged_payload_bytes"] == 40
    assert resumed.summary()["verified_received_payload_bytes"] == 20
    with pytest.raises(ValueError, match="already settled"):
        resumed.settle(retry, 15, status="complete")
    with pytest.raises(ValueError, match="settlement"):
        resumed.settle(interrupted, 21, status="complete")
    resumed.close()


def test_prior_receipts_are_explicit_identity_checked_and_counted_once(tmp_path):
    identity = {"source": "pinned"}
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps(dict(identity=identity, payload_bytes=90, range_requests=3)))
    ledger = TransferLedger(tmp_path / "ledger.sqlite", identity, prior_receipts=[receipt])
    assert ledger.summary()["charged_payload_bytes"] == 90
    ledger.close()
    resumed = TransferLedger(tmp_path / "ledger.sqlite", identity, prior_receipts=[receipt])
    assert resumed.summary()["charged_payload_bytes"] == 90
    resumed.close()
    with pytest.raises(ValueError, match="another source"):
        TransferLedger(tmp_path / "ledger.sqlite", {"source": "changed"})
    with pytest.raises(ValueError, match="Duplicate"):
        TransferLedger(tmp_path / "duplicate.sqlite", identity, prior_receipts=[receipt, receipt])
    with pytest.raises(ValueError, match="another source"):
        TransferLedger(tmp_path / "wrong.sqlite", {"source": "changed"}, prior_receipts=[receipt])


def test_scratch_cap_includes_replacement_temporary_and_sqlite_journal_capacity(tmp_path):
    scratch = ScratchBudget(tmp_path, 20000)
    output = tmp_path / "output.json"
    scratch.atomic_write(output, b"a" * 8192)
    with pytest.raises(LimitExceeded, match="disk limit"):
        scratch.atomic_write(output, b"b" * 12288)
    assert output.read_bytes() == b"a" * 8192
    assert not list(tmp_path.glob("*.tmp"))
    scratch.atomic_write(output, b"b" * 8192)
    assert scratch.receipt()["peak_observed_bytes"] >= 16384
    with pytest.raises(LimitExceeded, match="disk limit"):
        scratch.reserve_sqlite(tmp_path / "ledger.sqlite", 4096)
    with pytest.raises(ValueError, match="outside"):
        scratch.atomic_write(tmp_path.parent / "outside.json", b"outside")
