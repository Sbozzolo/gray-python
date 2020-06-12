![GRay-python-logo](gray-python.png)

# GRay-python (pre and post-processing)

![Tests](https://github.com/Sbozzolo/gray-python/workflows/Tests/badge.svg)
[![codecov](https://codecov.io/gh/Sbozzolo/gray-python/branch/master/graph/badge.svg?token=P8JZQHJWPC)](https://codecov.io/gh/Sbozzolo/gray-python)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

[GRay2](https://github.com/luxsrc/gray) is a general-relativistic
ray-tracing/radiation-transfer code based on the
[Lux](https://github.com/luxsrc/lux) framework. GRay2 can run anywhere: CPUs,
GPUs, Accelerators, and even FPGAs. With OpenCL, GRay2 always makes the best use
of the underlying hardware architecture, which results in the unparalleled
performance.

``GRay-python`` is a set of companion tools for Gray2. This package provides
modules to prepare GRay2 simulations and to analyze them. The documentation
for ``GRay-Python`` can be found 

## Tests

``GRay-python`` comes with a suite of unit tests. To run the tests, 
```sh
python3 -m unittest
```
Tests are automatically run after each commit by GitHub Actions.

## Documentation

``GRay-python`` uses Sphinx to generate the documentation. The extension
`syphinx-argparse` is needed. To produce the documentation
```sh
cd docs && make html
```
Documentation is automatically generated after each commit by GitHub Actions.
