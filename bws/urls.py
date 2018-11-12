from django.conf.urls import url
from rest_framework_swagger.views import SwaggerResourcesView, SwaggerApiView
from bws.view import BwsSwaggerUIView
from bws import vcf2prs, evaluation_timings, rest_api


urlpatterns = [
    url(r'^$', BwsSwaggerUIView.as_view(), name="django.swagger.base.view"),
    url(r'^api-docs/$', SwaggerResourcesView.as_view(), name="django.swagger.resources.view"),
    url(r'^api-docs/(?P<path>.*)/?$', SwaggerApiView.as_view(), name='django.swagger.api.view'),
]

internal_url_rest_patterns = [
    # url(r'^risk-factors/', risk_factors_ws.RiskFactorsView.as_view(), name='risk_factors'),
    url(r'^vcf2prs/', vcf2prs.Vcf2PrsView.as_view(), name='prs'),
    url(r'^evaluation/', evaluation_timings.EvaluationView.as_view(), name='evaluation'),
    url(r'^combine/', rest_api.CombineModelResultsView.as_view(), name='combine'),
]
