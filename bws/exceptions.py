from rest_framework.exceptions import APIException
from rest_framework import status
from django.utils.translation import ugettext_lazy as _
from rest_framework.exceptions import ValidationError


class PersonError(ValidationError):
    """
    PersonError raised for errors in the individuals input.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Person Error': detail})


class CancerError(ValidationError):
    """
    CancerError raised for errors in the cancer diagnosis input.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Cancer Error': detail})


class PathologyError(ValidationError):
    """
    PathologyError raised for errors in the pathology.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Pathology Error': detail})


class PedigreeFileError(ValidationError):
    """
    PedigreeFileError raised for errors in the input file.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Pedigree File Error': detail})


class PedigreeError(ValidationError):
    """
    PedigreeError raised for errors in the pedigree.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Pedigree Error': detail})


class GeneticTestError(ValidationError):
    """
    GeneticTestError raised for an error in a gene test.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Gene Test Error': detail})


class RiskFactorError(ValidationError):
    """
    RiskFactorError raised for an error in a gene test.
    @param detail: explanation of the error
    """

    def __init__(self, detail):
        super().__init__({'Risk Factor Error': detail})


class TimeOutException(APIException):
    ''' Request timeout '''
    status_code = status.HTTP_408_REQUEST_TIMEOUT
    default_detail = _('Request has timed out.')
    default_code = 'timeout'
