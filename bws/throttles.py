import logging
from rest_framework.throttling import UserRateThrottle, SimpleRateThrottle
from rest_framework.compat import is_authenticated

logger = logging.getLogger(__name__)


class LogThrottleMixin():
    """ Add logging for occasions when throttle stops request. """

    def allow_request(self, request, view):
        """ Implement the check to see if the request should be throttled. """
        success = super(LogThrottleMixin, self).allow_request(request, view)
        if not success:
            return self.throttle_fail(request.user)
        return success

    def throttle_fail(self, user):
        """ Called when a request to the API has failed due to throttling. """
        logger.warning("Request throttled (" + self.__class__.__name__ + "); USER: "+str(user))
        return False


class BurstRateThrottle(LogThrottleMixin, UserRateThrottle):
    """ Throttle short burst of requests from a user. """
    scope = 'burst'


class SustainedRateThrottle(LogThrottleMixin, UserRateThrottle):
    """ Throttle sustained requests from a user. """
    scope = 'sustained'


class EndUserIDRateThrottle(LogThrottleMixin, SimpleRateThrottle):
    """
    Limits the rate of API calls that may be made by a given end user.
    The end user id plus the user will be used as a unique cache key.
    """
    scope = 'enduser_burst'

    def get_cache_key(self, request, view):
        if is_authenticated(request.user):
            ident = str(request.user.pk)
        else:
            ident = str(self.get_ident(request))

        try:
            ident = ident+"$"+request.data.get('user_id')
            return self.cache_format % {
                'scope': self.scope,
                'ident': ident
            }
        except TypeError:
            return None
