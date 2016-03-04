======
Splipy
======

This repository contains the Splipy packages.


Dependencies
------------

This library depends on numpy and scipy. E.g. on Ubuntu::

    sudo apt-get install python-numpy python-scipy

If you use Python 2, you need OpenCV for the imaging tools to work::

    sudo apt-get install opencv


Installing
----------

To install, use::

    python setup.py install


To generate a package, use::

    python setup.py sdist --dist-dir .


Documenation
------------

To generate the documentation, first install sphinx. Then, in the doc folder::

    make html


Tests
-----

To run the tests, you can use your favourite test runner. For example, with
pytest::

    py.test splipy test_utils
