from rest_framework.exceptions import APIException
from rest_framework import status
from django.utils.translation import ugettext_lazy as _
from rest_framework.exceptions import ValidationError


class CanRiskError(ValidationError):

    def __init__(self, detail, famid=None):
        detail = ('['+famid+'] - ' if famid is not None and not "XXXX" else '') + detail
        super().__init__({self.__class__.err: detail})


class PersonError(CanRiskError):
    """
    PersonError raised for errors in the individuals input.
    @param detail: explanation of the error
    """
    err = 'Person Error'


class CancerError(CanRiskError):
    """
    CancerError raised for errors in the cancer diagnosis input.
    @param detail: explanation of the error
    """
    err = 'Cancer Error'


class PathologyError(CanRiskError):
    """
    PathologyError raised for errors in the pathology.
    @param detail: explanation of the error
    """
    err = 'Pathology Error'


class PedigreeFileError(CanRiskError):
    """
    PedigreeFileError raised for errors in the input file.
    @param detail: explanation of the error
    """
    err = 'Pedigree File Error'


class PedigreeError(CanRiskError):
    """
    PedigreeError raised for errors in the pedigree.
    @param detail: explanation of the error
    """
    err = 'Pedigree Error'


class GeneticTestError(CanRiskError):
    """
    GeneticTestError raised for an error in a gene test.
    @param detail: explanation of the error
    """
    err = 'Gene Test Error'


class RiskFactorError(CanRiskError):
    """
    RiskFactorError raised for an error in a gene test.
    @param detail: explanation of the error
    """
    err = 'Risk Factor Error'


class ModelError(ValidationError):
    """
    ModelError raised for an error when Fortran model code is run with non-zero exit code.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Model Error': detail})


class TimeOutException(APIException):
    ''' Request timeout '''
    status_code = status.HTTP_408_REQUEST_TIMEOUT
    default_detail = _('Request has timed out.')
    default_code = 'timeout'
