from collections import Counter
from typing import Optional


def common_caller(expt_callers: list[str]) -> Optional[str]:
    """
    Checks what variant caller was used across all experiments,
    where `expt_callers` is a list of used variant callers

    Raises a ValueError if more than one variant caller was used across experiments.
    Returns the common variant caller if all experiments used the same caller, or None if no callers
    """
    if len(expt_callers) == 0:
        return None
    caller_counts = Counter([caller for caller in expt_callers])
    if len(caller_counts) > 1:
        raise ValueError(
            "Found more than one variant caller used across experiments: "
            + f"{', '.join([f'{v} experiment(s) used {c}' for c, v in caller_counts.items()])}."
        )
    return caller_counts.most_common()[0][0]
