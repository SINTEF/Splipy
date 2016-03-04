======
Splipy
======

This repository contains the Splipy packages.


Dependencies
------------

This library depends on numpy and scipy. E.g. on Ubuntu::

    pip install numpy
    pip install scipy

If you use Python 2, you need OpenCV for the imaging tools to work::

    sudo apt-get install opencv

To generate the documentation you will need Sphinx::

    pip install sphinx

And to run the tests you can use your favourite test runner, for example
pytest::

    pip install pytest


Installing
----------

To install, use::

    python setup.py install


To generate a package, use::

    python setup.py sdist --dist-dir .


Documentation
-------------

To generate the documentation, run in the `doc` folder::

    make html


Tests
-----

To run the tests, you can use your favourite test runner. For example, with
pytest::

    py.test splipy test_utils
