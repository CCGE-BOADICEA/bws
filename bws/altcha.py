# ALTCHA server endpoints:
# GET /altcha - use this endpoint as challengeurl for the widget
# POST /submit - use this endpoint as the form action
# POST /submit_spam_filter - use this endpoint for form submissions with spam filtering
import logging
import os
import time

from altcha.altcha import ChallengeOptions, create_challenge, \
    verify_server_signature, verify_fields_hash, verify_solution
from django.conf import settings
from rest_framework.response import Response
from rest_framework.views import APIView
from drf_spectacular.utils import extend_schema


logger = logging.getLogger(__name__)


# Get HMAC key from environment variables
ALTCHA_HMAC_KEY = os.getenv(
    "ALTCHA_HMAC_KEY", settings.ALTCHA_HMAC_KEY
)


class ChallengeView(APIView):

    @extend_schema(exclude=True)    # exclude from the swagger docs
    def get(self, _request):
        try:
            challenge = create_challenge(
                ChallengeOptions(
                    hmac_key=ALTCHA_HMAC_KEY,
                    max_number=50000,
                )
            )
            return Response(challenge.to_dict())
        except Exception as e:
            return Response({"error": f"Failed to create challenge: {str(e)}"}, status=500)


class SubmitView(APIView):

    @extend_schema(exclude=True)    # exclude from the swagger docs
    def post(self, request):
        form_data = request.form.to_dict()
        payload = request.form.get("altcha")
        if not payload:
            return Response({"error": "Altcha payload missing"}, status=400)

        try:
            # Verify the solution
            verified, _err = verify_solution(payload, ALTCHA_HMAC_KEY, True)
            if not verified:
                return Response({"error": "Invalid Altcha payload"}, status=400)

            return Response({"success": True, "data": form_data})
        except Exception as e:
            return Response({"error": f"Failed to process Altcha payload: {str(e)}"}, status=400)


class SubmitSpamFilter(APIView):

    @extend_schema(exclude=True)    # exclude from the swagger docs
    def post(self, request):
        form_data = request.form.to_dict()
        payload = request.form.get("altcha")
        if not payload:
            return Response({"error": "Altcha payload missing"}, status=400)

        try:
            verified, verification_data, _err = verify_server_signature(
                payload, ALTCHA_HMAC_KEY
            )
            if not verified:
                return Response({"error": "Invalid Altcha payload"}, status=400)

            if verification_data.verified and int(verification_data.expire) > int(time.time()):
                if verification_data.classification == "BAD":
                    return Response({"error": "Classified as spam"}, status=400)

                if verification_data.fieldsHash:
                    verified = verify_fields_hash(
                        form_data,
                        verification_data.fields,
                        verification_data.fieldsHash,
                        "SHA-256",
                    )
                    if not verified:
                        return Response({"error": "Invalid fields hash"}, status=400)

                return Response(
                    {
                        "success": True,
                        "data": form_data,
                        "verificationData": verification_data.__dict__,
                    }
                )
            else:
                return Response({"error": "Invalid Altcha payload"}, status=400)
        except Exception as e:
            return Response({"error": f"Failed to process Altcha payload: {str(e)}"}, status=400)
