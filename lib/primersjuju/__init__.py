__version__ = "0.1.0"

class PrimersJuJuError(Exception):
    """general errors"""
    pass

class PrimersJuJuDataError(PrimersJuJuError):
    """error raised on bad input data"""
    pass
