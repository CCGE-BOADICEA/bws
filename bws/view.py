# from rest_framework_swagger.views import SwaggerUIView, get_full_base_path
# import rest_framework_swagger as rfs
# from django.utils.safestring import mark_safe
# from django.conf import settings
# import json
# from django.shortcuts import render_to_response
#
#
# class BwsSwaggerUIView(SwaggerUIView):
#     '''
#     Override get method to work with Django 1.11 which requires dict rather than RequestContext
#     docs.djangoproject.com/en/1.11/releases/1.11/#django-template-backends-django-template-render-prohibits-non-dict-context
#     '''
#
#     def get(self, request, *args, **kwargs):
#
#         if not self.has_permission(request):
#             return self.handle_permission_denied(request)
#
#         template_name = rfs.SWAGGER_SETTINGS.get('template_path')
#         data = {
#             'swagger_settings': {
#                 'discovery_url': "%s/api-docs/" % get_full_base_path(request),
#                 'api_key': rfs.SWAGGER_SETTINGS.get('api_key', ''),
#                 'api_version': rfs.SWAGGER_SETTINGS.get('api_version', ''),
#                 'token_type': rfs.SWAGGER_SETTINGS.get('token_type'),
#                 'enabled_methods': mark_safe(
#                     json.dumps(rfs.SWAGGER_SETTINGS.get('enabled_methods'))),
#                 'doc_expansion': rfs.SWAGGER_SETTINGS.get('doc_expansion', ''),
#             },
#             'rest_framework_settings': {
#                 'DEFAULT_VERSIONING_CLASS':
#                     settings.REST_FRAMEWORK.get('DEFAULT_VERSIONING_CLASS', '')
#                     if hasattr(settings, 'REST_FRAMEWORK') else None,
#
#             },
#             'django_settings': {
#                 'CSRF_COOKIE_NAME': mark_safe(
#                     json.dumps(getattr(settings, 'CSRF_COOKIE_NAME', 'csrftoken'))),
#             }
#         }
#         response = render_to_response(
#             template_name, data)
#
#         return response
