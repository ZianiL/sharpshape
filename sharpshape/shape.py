"""
Shape archetype of a pocket from its eigenvalue ratios rEV21 and rEV32.

The thresholds are the medians of rEV21 and rEV32 in the reference subset of the RNA pocket atlas
(278 pockets), used with all their decimals.
"""

MEDIAN_REV21 = 0.5896164508277054
MEDIAN_REV32 = 0.5432751373468785


def classify_shape(rev21, rev32):
    """Return 'sphere', 'disk', 'rod' or 'anisotropic' (strongly anisotropic)."""
    if rev21 >= MEDIAN_REV21:
        return "sphere" if rev32 >= MEDIAN_REV32 else "disk"
    return "rod" if rev32 >= MEDIAN_REV32 else "anisotropic"
