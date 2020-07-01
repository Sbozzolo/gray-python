Geodesic integration in a given spacetime
=========================================

``GRay-python`` has a module, ``gray.generate_metric``, to produce OpenCL
kernels that GRay2 can use to perform geodesic integration in any given
spacetime. This page describes how to use this functionality.

generate_gray_metric
------------------------

``GRay-python`` comes with a command-line tool ``generate_gray_metric``. This
executable takes as input a Python file and returns an OpenCL kernel ready to be
used in GRay2.

Input file and metric() function
________________________________

The input file has to be a valid Python3 code that defines a function
``metric()``. This function must take no arguments and must return a tuple with
the elements of the metric. The return value will probably look like:

.. code-block:: python

    return (g_tt, g_tx, g_ty, g_tz,
            g_xt, g_xx, g_xy, g_xz,
            g_yt, g_yx, g_yy, g_yz,
            g_zt, g_zx, g_zy, g_zz)

The most simple example is for Minkowski:

.. code-block:: python

    def metric():
        return (-1, 0, 0, 0,
                 0, 1, 0, 0
                 0, 0, 1, 0
                 0, 0, 0, 1)

A file containing only these lines would be a valid input for
``generate_gray_metric``.

If, like in most cases, the metric depends on the spacetime event, then the
input file must import ``sympy`` and must define the symbols that it uses.
For example,

.. code-block:: python

    from sympy import symbols
    t, x, y, z = symbols('t, x, y, z', real=True)

    def metric():
        return (-t + x + y, 0, 0, 0,
                0, x + y, 0, 0,
                0, 0, y, 0,
                0, 0, 0, z*z)


Note, the input file can define any helper function it may need, and can
introduce additional symbols. For instance:


.. code-block:: python

    from sympy import symbols
    t, x, y, z = symbols('t, x, y, z', real=True)
    a = symbols('a', real=True)

    def za():
        return z*a

    def metric():
        return (-t + x + y, 0, 0, 0,
                0, x + y, 0, 0,
                0, 0, y, 0,
                0, 0, 0, za())


``generate_gray_metric`` is turn all the additional symbols in run-time
adjustable arguments.

The concept of *levels*
________________________

The geodesic equation can be written in multiple different ways. The *best* way
to write the equation depends on the metric and on the hardware. For instance,
in some cases the analytical expression is more expensive to compute than the
numerical one. Or, when the metric is simple, it is possible to reduce
symbolically the expressions. In ``generate_gray_metric``, we implement many
different possible choices that we call *levels*.

Let us start from the levels 7, 8 and 9. In these levels, we deal directly with
the Christoffel symbols.

The most simple form of the geodesic equation is

.. math::
     \ddot{x}^\mu = -\Gamma^{\mu}_{\;\;\alpha\beta} \dot{x}^\alpha \dot{x}^\beta\,.
     :label: geo-eq-chri

Now, we have three choices. The first, is to compute symbolically the entire
right-hand-side of :eq:`geo-eq-chri`. In this case, ``generate_gray_metric``
computes analytically the Christoffel symbols, and use them to compute the
right-hand-side. Then, it turns the resulting expression into code. This is
level 7. A second possibility is to still analytically the Christoffel symbols,
but let GRay2 to compute numerically the dot products. This possibility is
implemented in level 8, and consists in using the language of matrices, so that
Equation :eq:`geo-eq-chri` reads

.. math::
  \ddot{x}^\mu = -\left(\vec{\Gamma}^{\mu} \vec{\dot{x}}\right) \cdot \vec{\dot{x}}\,.
  :label: geo-eq-chri-max

Finally, we can compute the derivatives analytically, and everything else
numerically (the Christoffel symbols and the dot products), as done in level 9.
Here, we use that the Christoffel symbols are

.. math::
  \Gamma^{\mu}_{\;\;\alpha\beta} = \frac{1}{2} g^{\mu\nu} \left(g_{\nu\alpha, \beta} + g_{\nu\beta, \alpha} - g_{\alpha\beta, \nu}\right)\,.
  :label: christoffel


The other levels are based on mathematical manipulations. With Equation
:eq:`christoffel`, Equation :eq:`geo-eq-chri` can be written as

.. math::
  \ddot{x}^\mu = - \left(g^{\mu\beta} \dot{x}^\alpha - \frac{1}{2}g^{\mu\alpha} \dot{x}^\beta  \right) g_{\beta\gamma, \alpha} \dot{x}^\gamma \,.
  :label: geo-eq-der

Defining

.. math::
   \begin{align}
   \Phi^\gamma_{\alpha \beta} &= g_{\beta\gamma, \alpha} \\
   \Upsilon &=  -\vec{\Phi}_{\gamma} \dot{x}^\gamma
   \end{align}
   :label: phi

Then,

.. math::
  \ddot{x}^\mu =  g^{\mu\beta}\Upsilon_{\alpha\beta} \dot{x}^\alpha + \frac{1}{2}g^{\mu\alpha}\Upsilon_{\alpha\beta} \dot{x}^\beta  \,.
  :label: geo-eq-der-ups

Or

.. math::
  \ddot{x}^\mu =  g^{\mu\beta}(\Upsilon^T)_{\beta\alpha} \dot{x}^\alpha + \frac{1}{2}g^{\mu\alpha}\Upsilon_{\alpha\beta} \dot{x}^\beta  \,.
  :label: geo-eq-der-ups-comp

In matrix notation

.. math::
  \ddot{\vec{x}} =  (\vec{g}^{-1}\vec{\Upsilon}^T) \cdot \dot{\vec{x}} + \frac{1}{2} (\vec{g}^{-1} \vec{\Upsilon}) \cdot \dot{\vec{x}}  \,.
  :label: geo-eq-der-ups-mat

Defining

.. math::
  \vec{\Xi} = \vec{\Upsilon}^T + \frac{1}{2} \vec{\Upsilon}
  :label: Xi

and

.. math::
  \vec{\Lambda} = \vec{g}^{-1} \cdot \vec{\Xi}
  :label: Lambda

The geodesic equation is simply

.. math::
  \ddot{\vec{x}} =  \vec{\Lambda} \dot{\vec{x}} \,.
  :label: geo-eq-Lambda

Levels 1-6 follow this derivation implementing part analytical and part
numerical computation. Increasing the level means that more and more
calculations are done numerically. Level 1 follows the entire derivation until
the end analytically, and provide a symbolical right-hand-side. Level 2 computes
:math:`\Lambda` analytically, and leave the final dot product to be computed
numerically. Level 3 is up to :math:`\Xi`, and the geodesic equation is computed numerically as 

.. math::
  \ddot{\vec{x}} = \vec{g}^{-1} \cdot \vec{\Xi} \cdot \dot{\vec{x}} \,.
  :label: geo-xi

At level 4, also :math:`\Xi` is computed numerically, with up :math:`\Upsilon`
analytical. Level 5, :math:`\Upsilon` is computed numerically from an analytical
:math:`\Phi`, and finally at Level 6, even :math:`\Phi` is computed numerically.


What level should I choose?
___________________________

If you have a very simple metric, choose levels with more analytical
computations. You could even add some ``simplify`` function calls in the code.
This, in general, increases the code generation time by several orders of
magnitude, and has benefits only if the metric is simple. In case you have
massively parallel devices, it is typically better to do more numerical
calculations.

Our suggestion is to try them all and find the fastest for the specific case in
consideration.



Command-line arguments
________________________

.. argparse::
   :filename: ../bin/generate_gray_metric
   :func: parser
   :prog: generate_gray_metric
   :nodescription:
