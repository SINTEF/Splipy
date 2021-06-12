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
---------
The library is packaged on splipy and can be installed through pip by simply running ::

    pip install splipy


Resources
---------

* `Getting-started guide <https://github.com/sintef/Splipy/tree/master/doc/Tutorial/Getting%20Started.ipynb>`_ - tutorial page (run with `jupyter <http://jupyter.org/>`_ to get interactive features)
* `Examples page <https://github.com/sintef/Splipy/tree/master/examples>`_ - stand-alone executables of different applications
* `API documentation <http://sintef.github.io/Splipy>`_ - technical details on all functions
* `Package installation page <https://pypi.org/project/Splipy>`_ - splipy on PyPi, packaged and ready for installation


======
Development and building from source
======

Dependencies
------------

**Required**

This library requires numpy and scipy. For building, cython is also
required. E.g. on Ubuntu::

    pip install numpy
    pip install scipy
    pip install cython

**Optional**

To use image processing tools, you need OpenCV ::

    pip install python-opencv

To generate the documentation you will need Sphinx::

    pip install sphinx

And to run the tests you can use your favourite test runner, for example
pytest::

    pip install pytest pytest-benchmark pytest-cov


Installing
----------

To install, use::

    python setup.py build_ext --inplace
    python setup.py install

To generate a package, use::

    python setup.py sdist --dist-dir .


Documentation
-------------

To generate the documentation, run in the `doc` folder::

    make html

To push generated docs online on the ``gh-pages`` branch, run the helper script::

    python push_documentation.py [remote]

where ``remote`` is the name of the remote to push to. If not given, it will be asked.


Tests
-----

To run the tests, you can use your favourite test runner. For example, with
pytest::

    PYTHONPATH=. py.test --benchmark-skip

To get a report of test coverage, run::

    PYTHONPATH=. py.test --benchmark-skip --cov=splipy --cov-report term-missing

Code analysis
-------------
You can use pylint3 to perform static code analysis on the module.
This can help identify bugs and give suggestions for improvements.

To install, use::

    pip3 install pylint

To perform the code analysis, use::

    pylint -d C --rcfile=pylint.cfg splipy/


Releasing
---------

To make a new release, it is recommended to install `bumpversion
<https://pypi.python.org/pypi/bumpversion>`_. To make a new release, run::

    bumpversion <type>

where `type` is one of `patch`, `minor` or `major`. This will up the version
number, create a commit and a tag. To push this to github, use::

    git push --tags

After that, Travis CI should automatically build and deploy the
packages to PyPi. It would be helpful to monitor the Travis build so
that errors can be fixed quickly. See the `list of builds
<https://travis-ci.org/sintefmath/Splipy/builds>`_.


======
Citations
======

If you use Splipy in your work, please consider citing
`K. A. Johannessen and E. Fonn 2020 J. Phys.: Conf. Ser. 1669 012032 <https://iopscience.iop.org/article/10.1088/1742-6596/1669/1/012032/meta>`_.
