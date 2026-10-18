"""SHAPE morphological archetype of a pocket: fixed median-threshold decision tree on rEV21 / rEV32.

Thresholds = medians of rEV21 and rEV32 over the non-redundant subset of the RNA pocket atlas (331 pockets),
used exactly (shown rounded as 0.5745 / 0.5433): the rounded values would move the median pocket to another archetype.
rEV31 = rEV21 * rEV32 is algebraically dependent and is not used. Not a trained classifier.
"""
MEDIAN_REV21 = 0.5744716108261702
MEDIAN_REV32 = 0.5433004248765029


def classify_shape(rev21, rev32):
    """'sphere', 'disc', 'rod' or 'anisotropic' (strongly anisotropic)."""
    if rev21 >= MEDIAN_REV21:
        return "sphere" if rev32 >= MEDIAN_REV32 else "disc"
    return "rod" if rev32 >= MEDIAN_REV32 else "anisotropic"
