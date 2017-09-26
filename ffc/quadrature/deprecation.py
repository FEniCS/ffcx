# -*- coding: utf-8 -*-

class QuadratureRepresentationDeprecationWarning(DeprecationWarning):
    """Warning about deprecation of quadrature representation"""
    def __init__(self):
        msg = """
*** ==================================================== ***
*** FFC: quadrature representaion is deprecated! It will ***
*** likely be removed in 2018.1.0 release. Use uflacs    ***
*** representation instead.                              ***
*** ==================================================== ***"""
        super(QuadratureRepresentationDeprecationWarning, self).__init__(msg)


# Be very annoying - never ignore the warning
import warnings
warnings.simplefilter('always', QuadratureRepresentationDeprecationWarning)
del warnings


def issue_deprecation_warning():
    """Issue warning about quadrature deprecation"""
    import warnings
    warnings.warn(QuadratureRepresentationDeprecationWarning(), stacklevel=2)
