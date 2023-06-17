import pytest
from nomadic.util.bed import BEDRecord
from nomadic.util.exceptions import BEDFormatError


# --------------------------------------------------------------------------------
# Tests for BEDRecord
#
# --------------------------------------------------------------------------------


def test_bedrecord_from_line_valid():
    line = "chr1\t100\t200\texample"
    record = BEDRecord.from_line(line)

    assert record.chrom == "chr1"
    assert record.start == 100
    assert record.end == 200
    assert record.name == "example"


def test_bedrecord_from_line_invalid_format():
    # Too few columns
    line = "chr1\t100\t200"
    with pytest.raises(BEDFormatError):
        BEDRecord.from_line(line)

    # Too many columns
    line = "chr1\t100\t200\texample\t0.45"
    with pytest.raises(BEDFormatError):
        BEDRecord.from_line(line)


def test_bedrecord_post_init_valid():
    record = BEDRecord("chr1", 100, 200, "example")
    assert record.start > 0
    assert record.end > 0
    assert record.start <= record.end


def test_bedrecord_post_init_invalid_start_end():
    with pytest.raises(BEDFormatError):
        BEDRecord("chr1", -100, 200, "example")

    with pytest.raises(BEDFormatError):
        BEDRecord("chr1", 200, 100, "example")
