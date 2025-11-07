"""
FDPS Animator - Animation tools for FDPS simulation output
"""

from .fdps_animator import (
    FDPSData,
    FDPSReader,
    FDPSTimeSeries,
    FDPSAnimator
)

__version__ = "0.1.0"
__all__ = ['FDPSData', 'FDPSReader', 'FDPSTimeSeries', 'FDPSAnimator']
