class MetadataFormatError(Exception):
    """Error in format or contents of a metadata file"""


class BEDFormatError(Exception):
    """Error in the format or contents of a BED file"""


class ReferenceGenomeMissingError(Exception):
    """Reference genome has not been downloaded"""
