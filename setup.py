from setuptools import setup

setup(
    name="GRay",
    version="0.1",
    author='Gabriele Bozzola',
    author_email='gabrielebozzola@arizona.edu',
    packages=['gray', 'gray.etc'],
    scripts=['bin/generate_gray_metric'],
    description='Pre and post-processing tools for GRay2.',
    install_requires=["numpy", "matplotlib", "sympy"],
    license='LICENSE.txt',
)
