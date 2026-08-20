# ALTCHA server endpoint:
# GET /altcha - use this endpoint as the challenge attribute for the widget
import datetime
import logging

from altcha import create_challenge
from django.conf import settings
from rest_framework.response import Response
from rest_framework.views import APIView
from drf_spectacular.utils import extend_schema


logger = logging.getLogger(__name__)


class ChallengeView(APIView):

    @extend_schema(exclude=True)    # exclude from the swagger docs
    def get(self, _request):
        '''
        Fetches a new random proof-of-work (v2) challenge to be used by the ALTCHA widget.
        The challenge is signed with the shared HMAC key so that the solution posted back
        with a form can be verified as one this site issued (see AltchaFormMixin).
        '''
        try:
            expires_at = (datetime.datetime.now(datetime.timezone.utc) +
                          datetime.timedelta(seconds=settings.ALTCHA_EXPIRY_SECONDS))
            challenge = create_challenge(
                algorithm=settings.ALTCHA_ALGORITHM,
                cost=settings.ALTCHA_COST,
                expires_at=expires_at,
                hmac_secret=settings.ALTCHA_HMAC_KEY,
            )
            return Response(challenge.to_dict())
        except Exception as e:
            logger.error("Failed to create ALTCHA challenge: %s", e)
            return Response({"error": f"Failed to create challenge: {str(e)}"}, status=500)
