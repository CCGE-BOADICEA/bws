====
bws
====


bws is a Django app to run web-services for BOADICEA.

Quick start
-----------

1. Installation::

pip install -e git+ssh://git@github.com/CCGE-BOADICEA/bws.git#egg=bws

2. Add "bws" to your ``INSTALLED_APPS`` in ``settings.py``::

    INSTALLED_APPS = (
        ...
        'bws',
    )

3. Import the bws settings in ``settings.py``::

    from bws.settings import *

