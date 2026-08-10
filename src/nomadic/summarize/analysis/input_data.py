from collections import Counter


def common_caller(expt_callers: list[str]) -> str:
    """
    Checks what variant caller was used across all experiments,
    where `expt_callers` is a list of used variant callers

    Raises a ValueError if more than one variant caller was used across experiments.
    Returns the common variant caller if all experiments used the same caller
    """
    if len(expt_callers) == 0:
        raise ValueError("No experiment callers provided")
    caller_counts = Counter([caller for caller in expt_callers])
    if len(caller_counts) > 1:
        raise ValueError(
            "Found more than one variant caller used across experiments: "
            + f"{', '.join([f'{v} experiment(s) used {c}' for c, v in caller_counts.items()])}."
        )
    return caller_counts.most_common()[0][0]


def common_reference_name(reference_names: list[str]) -> str:
    """
    Checks what reference was used across all experiments,
    where `reference_names` is a list of used references

    Raises a ValueError if more than one reference was used across experiments.
    Returns the common reference name if all experiments used the same reference
    """
    if len(reference_names) == 0:
        raise ValueError("No experiment reference names provided")
    reference_counts = Counter([ref for ref in reference_names])
    if len(reference_counts) > 1:
        raise ValueError(
            "Found more than one reference used across experiments: "
            + f"{', '.join([f'{v} experiment(s) used {c}' for c, v in reference_counts.items()])}."
        )
    return reference_counts.most_common()[0][0]
