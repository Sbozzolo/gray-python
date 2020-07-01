Overview
========

GRay2_ is a general-relativistic ray-tracing/radiation-transfer code based on
the Lux_ framework. GRay2 can run anywhere: CPUs, GPUs, Accelerators, and even
FPGAs. With OpenCL, GRay2 always makes the best use of the hardware
architecture, which results in the unparalleled performance.

``GRay-python`` is a set of companion tools for Gray2. This package provides
modules to prepare GRay2 simulations and to analyze them. The documentation in
these pages concerns only ``GRay-Python``. For help with GRay2, see its own
documentation.

.. _GRay2: https://github.com/luxsrc/gray
.. _Lux: https://github.com/luxsrc/lux


Features
--------

Features currently implemented:

- Generate GRay2 kernels for arbitrary analytical spacetimes (``gray.generate_metric``)

Installation
------------

Clone the repo_:

.. _repo: https://github.com/Sbozzolo/gray-python:

.. code-block:: bash

   git clone https://github.com/Sbozzolo/gray-python.git

Move into the folder and install with pip:

.. code-block:: bash

   cd gray-python && pip3 install --user .

Usage
-----

.. toctree::
   :maxdepth: 1

   generate_metric.rst
