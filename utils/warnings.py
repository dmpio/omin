"""Warning tools."""

import warnings


def deprecated(func):
    """A decorator which can be used to mark functions as deprecated.

    It will result in a warning being emmitted when the function is used.

    """
    def newFunc(*args, **kwargs):
        warnings.warn("".join(["Call to deprecated function ", func.__name__]),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc
