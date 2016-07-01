''' Pagination, filter and elastic list and retrieve views. '''

from rest_framework.filters import DjangoFilterBackend, OrderingFilter
from rest_framework.response import Response
from rest_framework.pagination import LimitOffsetPagination
from django.http.response import Http404
from bws.elastic_obj import ElasticObject
import datetime


# class ElasticLimitOffsetPagination(LimitOffsetPagination):
#     ''' Extend L{LimitOffsetPagination} for pagination of elastic resources. '''
# 
#     def paginate_queryset(self, queryset, request, view=None):
#         if not hasattr(view, 'es_count'):
#             return super().paginate_queryset(queryset, request, view)
# 
#         self.limit = self.get_limit(request)
#         if self.limit is None:
#             return None
# 
#         self.offset = self.get_offset(request)
#         self.count = view.es_count
#         self.request = request
#         if self.count > self.limit and self.template is not None:
#             self.display_page_controls = True
#         return queryset
# 
#     def get_limit(self, request):
#         ''' Override so that if a limit per page is not set a default of 10 is used. '''
#         limit = super().get_limit(request)
#         if limit is None:
#             return 10
#         return limit


class BwsFilterBackend(OrderingFilter, DjangoFilterBackend):
    ''' Extend L{DjangoFilterBackend} for filtering elastic resources. '''

    def filter_queryset(self, request, queryset, view):
        ''' Override this method to request just the documents required from elastic. '''
#         q_size = view.paginator.get_limit(request)
#         q_from = view.paginator.get_offset(request)

        filterable = getattr(view, 'filter_fields', [])
        print(filterable)
        filters = dict([(k, v) for k, v in request.GET.items() if k in filterable])
        print(filters)
        mut_freq = filters.get('mut_freq', 'UK')
        brca1_mut_search_sensitivity = filters.get('brca1_mut_search_sensitivity', .7)
        brca2_mut_search_sensitivity = filters.get('brca1_mut_search_sensitivity', .8)

        print("XXXXXXXXXXXXXXXX1")
        print(request.data)
        print("XXXXXXXXXXXXXXXX2")

        results = []
        new_obj = {
            'mut_prob': {
                'brca1': 2.1,
                'brca2': 2.3,
                'no_mutation': 95.5
            },
            'date': datetime.datetime.now(),
            'model_params': {
                'family member': 'PB1',
                'mutation frequencies': 'UK, BRCA1: 6.394d-4, BRCA2: 0.00102',
                'mutation sensitivities': 'Default, BRCA1: 0.7, BRCA2: 0.8',
                'cancer incidence rated': 'UK'
            }
        }
        results.append(new_obj)

        new_obj = {
            'mut_prob': {
                'brca1': 2.0,
                'brca2': 2.0,
                'no_mutation': 96.0
            },
            'date': datetime.datetime.now(),
            'model_params': {
                'family member': 'PB1',
                'mutation frequencies': 'UK, BRCA1: 6.394d-4, BRCA2: 0.00102',
                'mutation sensitivities': 'Default, BRCA1: 0.7, BRCA2: 0.8',
                'cancer incidence rated': 'UK'
            }
        }
        results.append(new_obj)

        return results


class ListBwsMixin(object):
    ''' List queryset. '''
    filter_backends = [BwsFilterBackend, ]

    def get_queryset(self):
        return None

    def list(self, request, *args, **kwargs):
        ''' Retrieve a list of documents. '''
        qs = self.filter_queryset(self.get_queryset())
        page = self.paginate_queryset(qs)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        serializer = self.get_serializer(qs, many=True)
        return Response(serializer.data)


# class RetrieveBwsMixin(object):
#     ''' Retrieve an instance. '''
#     def retrieve(self, request, *args, **kwargs):
#         ''' Retrieve a document instance. '''
#         instance = self.get_object()
#         serializer = self.get_serializer(instance)
#         return Response(serializer.data)
# 
#     def get_object(self):
#         new_obj = ElasticObject(initial={
#                         'mut_prob': {
#                             'brca1': 2.1,
#                             'brca2': 2.3,
#                             'no_mutation': 95.5
#                         },
#                         'date': datetime.datetime.now(),
#                         'model_params': {
#                             'family member': 'PB1',
#                             'mutation frequencies': 'UK, BRCA1: 6.394d-4, BRCA2: 0.00102',
#                             'mutation sensitivities': 'Default, BRCA1: 0.7, BRCA2: 0.8',
#                             'cancer incidence rated': 'UK'
#                         }
#                     })
#         new_obj.uuid = 'fam_id1'
#         return new_obj
