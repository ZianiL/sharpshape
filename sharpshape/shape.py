"""Shape archetype of a pocket: fixed median-threshold decision tree on rEV21 / rEV32.

Thresholds = medians of rEV21 and rEV32 over the reference subset of the RNA pocket atlas (278 pockets),
used exactly (shown rounded as 0.5896 / 0.5433): rounded values could move pockets close to a threshold.
rEV31 = rEV21 * rEV32 is algebraically dependent and is not used. Not a trained classifier.
"""
MEDIAN_REV21 = 0.5896164508277054
MEDIAN_REV32 = 0.5432751373468785


def classify_shape(rev21, rev32):
    """'sphere', 'disk', 'rod' or 'anisotropic' (strongly anisotropic)."""
    if rev21 >= MEDIAN_REV21:
        return "sphere" if rev32 >= MEDIAN_REV32 else "disk"
    return "rod" if rev32 >= MEDIAN_REV32 else "anisotropic"
