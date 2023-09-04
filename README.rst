====
bws
====


BWS is a Django app to run web-services for `BOADICEA <https://canrisk.org/about/>`_.
It use the Django REST framework to provide the Web APIs.

Quick start
-----------

1. Installation::

    pip install -e git+https://git@github.com/CCGE-BOADICEA/bws.git#egg=bws
    

This will install the bws app and the dependencies.    

2. Locate the installation directory (e.g. with pip show bws). Update the ``bws.settings.py``.
In particular change the ``FORTRAN_HOME`` parameter to set the directory location for the cancer risk models.
Depending on the file structure it may be necessary to also change the ``HOME`` location in ``BC_MODEL``
and ``OC_MODEL`` in this file. These and the ``PROBS_EXE`` and ``RISKS_EXE`` parameters define the location
of the executables for the mutation probability and risk calculation::

    FORTRAN_HOME = "/usr/src/"
    ....
    
    #
    # BREAST CANCER MODEL
    BC_MODEL = {
        'NAME': 'BC',
        'HOME': os.path.join(FORTRAN_HOME, 'boadicea'),
        'PROBS_EXE': 'boadicea_probs.exe',
        'RISKS_EXE': 'boadicea_risks.exe',
    ....

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

    from bws import rest_api
    from rest_framework.authtoken.views import ObtainAuthToken
    ....
     
    url_rest_patterns = [
        path('boadicea/', rest_api.BwsView.as_view(), name='bws'),    # breast cancer risk model
        path('ovarian/', rest_api.OwsView.as_view(), name='ows'),     # ovarian cancer risk model
        path('auth-token/', ObtainAuthToken.as_view()),
        path('combine/', rest_api.CombineModelResultsView.as_view(), name='combine'),
    ]
    urlpatterns.extend(url_rest_patterns)


7. Run tests::

    python manage.py test bws.tests.test_bws \
                          bws.tests.test_ows.OwsTests \
                          bws.tests.test_throttling \
                          bws.tests.tests_pedigree_validation

The `run_webservice.py <https://github.com/CCGE-BOADICEA/bws/blob/master/bws/scripts/run_webservice.py>`_ 
script can be used to submit requests to the BWS web-service. It takes a username and
pedigree file and will prompt for the password and run the risk calculation via the web-service::

    ${PATH_TO_BWS}/bws/scripts/run_webservice.py --help
    ${PATH_TO_BWS}/bws/scripts/run_webservice.py --url ${URL} -u ${USER} \
                                                 -p ${PATH_TO_BWS}/bws/tests/data/pedigree_data.txt

License
-------

This software is distributed under the GPL either version 3 or later license.

Publication
-----------

pedigreejs: a web-based graphical pedigree editor. Carver T, et al. `Bioinformatics, Volume 34, Issue 6, 15 March 2018 <http://dx.doi.org/10.1093/bioinformatics/btx705>`_.
