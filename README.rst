====
bws
====


bws is a Django app to run web-services for BOADICEA.

Quick start
-----------

1. Installation::

    pip install -e git+ssh://git@github.com/CCGE-BOADICEA/bws.git#egg=bws

2. Add "rest_framework", "rest_framework.authtoken" and "bws" to your ``INSTALLED_APPS`` in ``settings.py``::

    INSTALLED_APPS = (
        ...
        'rest_framework',
        'rest_framework.authtoken',
        'bws',
    )

3. Import the bws settings in ``settings.py``::

    from bws.settings import *

4. Add settings for web-services throttling in ``settings.py``::

    REST_FRAMEWORK = {
        'DEFAULT_THROTTLE_RATES': {
            'sustained': '5000/day',
            'burst': '60/min',
            'enduser_burst': '50/min'
        }
    }

5. Add a working directory in ``settings.py``::

    CWD_DIR = "/tmp"

6. Run tests::

    python manage.py test bws.tests.test_bws \
                          bws.tests.test_ows.OwsTests \
                          bws.tests.test_throttling \
                          bws.tests.tests_pedigree_validation
