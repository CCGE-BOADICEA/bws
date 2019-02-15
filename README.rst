====
bws
====


bws is a Django app to run web-services for `BOADICEA <https://canrisk.org/about/>`_.

Quick start
-----------

1. Installation::

    pip install -e git+ssh://git@github.com/CCGE-BOADICEA/bws.git#egg=bws

2. Update the ``bws.settings.py``. In particular change the ``FORTRAN_HOME`` parameter to
set the directory location for the cancer risk models. Depending on the file structure
it may be necessary to also change the ``HOME`` location in ``BC_MODEL`` and ``OC_MODEL``
in this file. These and the ``PROBS_EXE`` and ``RISKS_EXE`` parameters define the location
of the executables for the mutation probability and risk calculation.

3. If you need to start a Django project::

    django-admin startproject [project_name]

4. Add "rest_framework", "rest_framework.authtoken" and "bws" to your ``INSTALLED_APPS`` in ``settings.py``::

    INSTALLED_APPS = (
        ...
        'rest_framework',
        'rest_framework.authtoken',
        'bws',
    )

5. Import the bws settings in ``settings.py``::

    from bws.settings import *
  
6. Add web-service endpoints to the ```urls.py``::

	 url_rest_patterns = [
	     url(r'^boadicea/', rest_api.BwsView.as_view(), name='bws'),
	     url(r'^ovarian/', rest_api.OwsView.as_view(), name='ows'),
	 ]
	 urlpatterns.extend(url_rest_patterns)

7. Run tests::

    python manage.py test bws.tests.test_bws \
                          bws.tests.test_ows.OwsTests \
                          bws.tests.test_throttling \
                          bws.tests.tests_pedigree_validation
