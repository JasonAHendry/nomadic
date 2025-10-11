import pytest
from nomadic.util.metadata import correct_barcode_format, MetadataTableParser
from nomadic.util.exceptions import MetadataFormatError


# --------------------------------------------------------------------------------
# Tests for check_barcode_format
#
# --------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "barcode,try_to_fix,expected",
    [
        ("barcode01", True, "barcode01"),  # Valid barcode
        ("barcode25", False, "barcode25"),  # Valid barcode, no fix needed
    ],
)
def test_check_barcode_format_good(barcode, try_to_fix, expected):
    assert correct_barcode_format(barcode, try_to_fix) == expected


@pytest.mark.parametrize(
    "barcode,try_to_fix,expected",
    [
        (26, True, "barcode26"),  # Fixable
        ("barcode_1", True, "barcode01"),  # Fixable
    ],
)
def test_check_barcode_warning(barcode, try_to_fix, expected):
    with pytest.warns(UserWarning):
        assert correct_barcode_format(barcode, try_to_fix) == expected


@pytest.mark.parametrize(
    "barcode,try_to_fix,expected",
    [
        (200, True, "barcode26"),
        ("barcode", True, "barcode01"),
        ("barcode01barcode02", True, "barcode01"),
    ],
)
def test_check_barcode_warning_error(barcode, try_to_fix, expected):
    with pytest.warns(UserWarning):
        with pytest.raises(MetadataFormatError):
            assert correct_barcode_format(barcode, try_to_fix) == expected


# --------------------------------------------------------------------------------
# Tests for MetadataTableParser
#
# --------------------------------------------------------------------------------

test_files_folder = "src/nomadic/util/_test_data/"


@pytest.mark.parametrize(
    "csv_path,csv_shape",
    [
        (test_files_folder + "metadata/sample_info_good.csv", (5, 5)),
        (test_files_folder + "metadata/sample_info_semicolon.csv", (5, 5)),
        (test_files_folder + "metadata/sample_info_eurosep.csv", (5, 6)),
        (test_files_folder + "metadata/sample_info_hybridsep-col.csv", (5, 6)),
    ],
)
def test_metadata_correct(csv_path, csv_shape):
    metadata = MetadataTableParser(csv_path)
    assert metadata.df.shape == csv_shape


@pytest.mark.parametrize(
    "csv_path,csv_shape",
    [
        (test_files_folder + "metadata/sample_info_badbarcode-format.csv", (5, 5)),
        (test_files_folder + "metadata/sample_info_badbarcode-int.csv", (5, 5)),
    ],
)
def test_metadata_warns(csv_path, csv_shape):
    with pytest.warns(UserWarning):
        metadata = MetadataTableParser(csv_path)
    assert metadata.df.shape == csv_shape


@pytest.mark.parametrize(
    "csv_path,error_msg",
    [
        (
            test_files_folder + "metadata/sample_info_onecolumn.csv",
            "Metadata must contain column called sample_id!",
        ),
        (
            test_files_folder + "metadata/sample_info_dupbarcode.csv",
            "Column barcode must contain only unique entires, but barcode04 is duplicated.",
        ),
        (
            test_files_folder + "metadata/sample_info_nobarcode.csv",
            "Metadata must contain column called barcode!",
        ),
        (
            test_files_folder + "metadata/sample_info_badheader.csv",
            "Found multiple delimiters (, ;) in header: ﻿barcode;sample_id,parasite_per_ul,location;platform;number.",
        ),
    ],
)
def test_metadata_errors(csv_path, error_msg):
    with pytest.raises(MetadataFormatError) as excinfo:
        _ = MetadataTableParser(csv_path)
    assert str(excinfo.value) == error_msg


@pytest.mark.parametrize(
    "csv_path",
    [
        test_files_folder + "metadata/sample_info_column_correction_case.csv",
        test_files_folder + "metadata/sample_info_column_correction_plural.csv",
        test_files_folder + "metadata/sample_info_column_correction_other.csv",
        test_files_folder + "metadata/sample_info_column_correction_space.csv",
        test_files_folder + "metadata/sample_info_column_correction_underscore.csv",
    ],
)
def test_metadata_column_corrections(csv_path):
    meta_table = MetadataTableParser(csv_path)
    assert "barcode" in meta_table.df.columns
    assert "sample_id" in meta_table.df.columns

    assert meta_table.df["barcode"].tolist() == [
        "barcode01",
        "barcode02",
        "barcode03",
        "barcode04",
        "barcode05",
    ]
    assert meta_table.df["sample_id"].tolist() == [
        "3D7 (NB01)",
        "Dd2 (NB02)",
        "CamC580Y (NB03)",
        "GB4 (NB04)",
        "HB3 (NB05)",
    ]
