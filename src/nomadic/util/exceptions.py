class NomadicError(Exception):
    """Base class for all Nomadic-related errors."""


class UserInputError(NomadicError):
    """Error due to invalid user input.

    These are errors that the user can correct by providing valid input.
    We catch those and show a user-friendly error message.
    """


class MetadataFormatError(UserInputError):
    """Error in format or contents of a metadata file"""


class BEDFormatError(UserInputError):
    """Error in the format or contents of a BED file"""


class ReferenceGenomeMissingError(UserInputError):
    """Reference genome has not been downloaded"""
