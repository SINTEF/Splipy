======
Splipy
======

This repository contains the Splipy packages. Splipy is a pure python library
for the creation, evaluation and manipulation of B-spline and NURBS geometries.
It supports n-variate splines of any dimension, but emphasis is made on the
use of curves, surfaces and volumes. The library is designed primarily for
analysis use, and therefore allows fine-grained control over many aspects which
is not possible to achieve with conventional CAD tools. 


Resources
---------

* `Getting-started guide <https://github.com/sintefmath/Splipy/tree/master/doc/Tutorial/Getting%20Started.ipynb>`_ - tutorial page (run with `jupyter <http://jupyter.org/>`_ to get interactive features)
* `Examples page <https://github.com/sintefmath/Splipy/tree/master/examples>`_ - stand-alone executables of different applications
* `API documentation <http://sintefmath.github.io/Splipy>`_ - technical details on all functions


Dependencies
------------

**Required**

This library requires numpy and scipy. E.g. on Ubuntu::

    pip install numpy
    pip install scipy

**Optional**

To use image processing tools, you need OpenCV ::

    pip install python-opencv

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

After that, to create the actual packages, run::

    rm -rf dist
    python setup.py sdist
    python setup.py bdist_wheel --universal

to create a source distribution and a wheel. These can then be uploaded where
they need to be uploaded. The recommended way to do that is using `twine
<https://pypi.python.org/pypi/twine>`_::

    twine upload dist/* -r <index>

Where `index` is the name of the index in your `~/.pypirc` where you want to
upload.
