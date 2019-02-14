import os
from setuptools import setup, find_packages

# with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as readme:
#     README = readme.read()

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))
ROOT = os.path.abspath(os.path.dirname(__file__))

setup(
    name='bws',
    version='0.0.1',
    packages=find_packages(),
    package_data={'bws': ['tests/data/*txt'], },
    include_package_data=True,
    zip_safe=False,
    url='https://github.com/CCGE-BOADICEA/bws',
    description='A Django app for web-services for BOADICEA.',
    long_description=open(os.path.join(ROOT, 'README.rst')).read(),
    install_requires=["requests>=2.7.0", "Django>=1.11,<2.0", "djangorestframework==3.9.0",
                      "coreapi==2.3.3"],
    classifiers=[
        'Environment :: Web Environment',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.4',
        'Topic :: Internet :: WWW/HTTP',
        'Topic :: Internet :: WWW/HTTP :: Dynamic Content',
    ],
)