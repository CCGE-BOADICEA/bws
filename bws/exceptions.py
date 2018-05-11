from rest_framework.exceptions import APIException
from rest_framework import status
from django.utils.translation import ugettext_lazy as _


class TimeOutException(APIException):
    ''' Request timeout '''
    status_code = status.HTTP_408_REQUEST_TIMEOUT
    default_detail = _('Request has timed out.')
    default_code = 'timeout'
