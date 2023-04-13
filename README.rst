.. image:: images/logo_small.svg

======
Splipy
======

This repository contains the Splipy packages. Splipy is a pure python library
for the creation, evaluation and manipulation of B-spline and NURBS geometries.
It supports n-variate splines of any dimension, but emphasis is made on the
use of curves, surfaces and volumes. The library is designed primarily for
analysis use, and therefore allows fine-grained control over many aspects which
is not possible to achieve with conventional CAD tools.


Installation
------------
The library is packaged on PyPI and can be installed through pip by simply
running ::

    pip install splipy


Resources
---------

* `Getting-started guide <https://github.com/sintef/Splipy/tree/master/doc/Tutorial/Getting%20Started.ipynb>`_ - tutorial page (run with `jupyter <http://jupyter.org/>`_ to get interactive features)
* `Examples page <https://github.com/sintef/Splipy/tree/master/examples>`_ - stand-alone executables of different applications
* `API documentation <http://sintef.github.io/Splipy>`_ - technical details on all functions
* `Package installation page <https://pypi.org/project/Splipy>`_ - splipy on PyPi, packaged and ready for installation


====================================
Development and building from source
====================================

Poetry
------

Splipy uses Poetry as a project management tool. To install poetry, use::

    pip install poetry

Poetry is the only tool that must be installed outside of the virtual
environment for Splipy. Once installed, run the command::

    make install

in the root Splipy directory. This will install Splipy and its dependencies in a
virtual environment located in the ``.venv`` directory.

To run the tests::

    make pytest


Installing
----------

To install, use::

    pip install .

To generate a package (source distribution or wheel), use::

    poetry build -f sdist
    poetry build -f wheel

Don't upload wheels to PyPI manually. They are built by CI runners whenever a
new version is tagged (see below).


Documentation
-------------

To generate the documentation, run::

    make doc

To push generated docs online on the ``gh-pages`` branch, run the helper script::

    python push_documentation.py [remote]

where ``remote`` is the name of the remote to push to. If not given, it will be asked.


Tests
-----

To run the tests, use::

    make pytest

For benchmarks::

    make bench


Code analysis
-------------

You can use pylint3 to perform static code analysis on the module.
This can help identify bugs and give suggestions for improvements.

To install, use::

    poetry run pip install pylint

To perform the code analysis, use::

    poetry run pylint -d C --rcfile=pylint.cfg splipy/


Releasing
---------

To make a new release, it is recommended to install `bumpversion
<https://pypi.python.org/pypi/bumpversion>`_. To make a new release, run::

    bumpversion <type>

where `type` is one of `patch`, `minor` or `major`. This will up the version
number, create a commit and a tag. To push this to github, use::

    git push --tags

After that, CI should automatically build and deploy the packages to pyPi. It
would be helpful to monitor the actions so that errors can be fixed quickly.


=========
Citations
=========

If you use Splipy in your work, please consider citing
`K. A. Johannessen and E. Fonn 2020 J. Phys.: Conf. Ser. 1669 012032 <https://iopscience.iop.org/article/10.1088/1742-6596/1669/1/012032/meta>`_.
