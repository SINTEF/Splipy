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

You should activate this virtual environment whenever you work on Splipy. The
makefile commands do not require it, but it's a good habit::

    source .venv/bin/activate

To run the tests::

    make test


Installing
----------

To install, use::

    pip install .

To generate a package (source distribution or wheel), use::

    make sdist
    make wheel
    make build  # both sdist and wheel

Don't upload wheels to PyPI manually. They are built by CI runners whenever a
new version is tagged (see below).


Tests
-----

To run the tests, use::

    make test

To run specific parts of the tests, use::

    make pytest
    make mypy
    make lint-check

The lint-check stage of the tests will complain about linter and style errors.
Some of these can be fixed automatically. To do this, run::

    make lint
    make format

For benchmarks::

    make benchmark


Documentation
-------------

To generate the documentation, run::

    make doc

To push generated docs online on the ``gh-pages`` branch, run the helper script::

    python push_documentation.py [remote]

where ``remote`` is the name of the remote to push to. If not given, it will be asked.


Releasing
---------

To make a new release, run the `bump-my-version` command::

    bump-my-version --dry-run <part>

Where `<part>` is the part you want to bump: either `major`, `minor`, `patch`,
`pre_label` or `pre_number`.

You can also specify the new version directly by using::

    bump-my-version --dry-run --new-version <new_version>

Once you are satisfied with the results, run the command without `--dry-run`.
We highly recommend to always use a dry run!

After that, CI should automatically build and deploy the packages to PyPi. It
would be helpful to monitor the actions so that errors can be fixed quickly.


=========
Citations
=========

If you use Splipy in your work, please consider citing
`K. A. Johannessen and E. Fonn 2020 J. Phys.: Conf. Ser. 1669 012032 <https://iopscience.iop.org/article/10.1088/1742-6596/1669/1/012032/meta>`_.
