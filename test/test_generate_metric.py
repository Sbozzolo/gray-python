#!/usr/bin/env python3

import unittest
import unittest.mock
import sympy
from gray import generate_metric
from importlib.resources import read_text  # Used to load test results


class TestGenerateMetric(unittest.TestCase):

    def setUp(self):
        self.t, self.x, self.y, self.z, self.a = sympy.symbols('t, x, y, z, a', real=True)
        self.a_spin = sympy.symbols('a_spin', real=True)
        self.uUP = sympy.symbols('uUPt, uUPx, uUPy, uUPz', real=True)

        self.M = sympy.Matrix([[-self.t * self.a, 0, 0, self.t],
                               [0, self.x, 0, 0],
                               [0, 0, self.y, 0],
                               [self.t, 0, 0, self.z **2]])

        self.g_deriv = sympy.derive_by_array(self.M, [self.t, self.x, self.y, self.z])

        # Kerr Schild
        aa = self.a_spin * self.a_spin
        zz = self.z * self.z
        kk = 1 / 2 * (self.x * self.x + self.y * self.y + zz - aa)
        dd = sympy.sqrt(kk * kk + aa * zz)
        rr = dd + kk
        r = sympy.sqrt(rr)

        f = 2 * rr * r / (rr * rr + aa * zz)
        lx = (r * self.x + self.a_spin * self.y) / (rr + aa)
        ly = (r * self.y - self.a_spin * self.x) / (rr + aa)
        lz = self.z / r

        g_tt, g_tx, g_ty, g_tz = -1 + f * 1 * 1, f * 1 * lx, f * 1 * ly, f * 1 * lz
        g_xt, g_xx, g_xy, g_xz = f * lx * 1, 1 + f * lx * lx, f * lx * ly, f * lx * lz
        g_yt, g_yx, g_yy, g_yz = f * ly * 1, f * ly * lx, 1 + f * ly * ly, f * ly * lz
        g_zt, g_zx, g_zy, g_zz = f * lz * 1, f * lz * lx, f * lz * ly, 1 + f * lz * lz

        self.g_ks = sympy.Matrix([[g_tt, g_tx, g_ty, g_tz], [g_xt, g_xx, g_xy, g_xz],
                                  [g_yt, g_yx, g_yy, g_yz,], [g_zt, g_zx, g_zy, g_zz]])

        self.g_ksUP = sympy.inv_quick(self.g_ks)
        self.g_ks_deriv = sympy.derive_by_array(self.g_ks,
                                                [self.t, self.x, self.y, self.z])

    def test_metric_from_tuple(self):

        def metric():
            return (-self.t * self.a, 0, 0, self.t,
                    0, self.x, 0, 0,
                    0, 0, self.y, 0,
                    0, 0, 0, self.z ** 2)

        with self.assertRaises(ValueError):
            generate_metric.metric_from_tuple((1, 2))

        with self.assertRaises(TypeError):
            generate_metric.metric_from_tuple(1)

        self.assertEqual(generate_metric.metric_from_tuple(metric()),
                         self.M)

    def test_find_indeces_for_contraction(self):

        self.assertEqual(generate_metric.find_indeces_for_contraction((0, 5), (2, 3)),
                         (1, 2))
        self.assertEqual(generate_metric.find_indeces_for_contraction((0, 1), (2, 3)),
                         (0, 1))
        self.assertEqual(generate_metric.find_indeces_for_contraction((1, 0), (2, 3)),
                         (0, 1))
        self.assertEqual(generate_metric.find_indeces_for_contraction((3, 2), (0, 1)),
                         (0, 1))
        self.assertEqual(generate_metric.find_indeces_for_contraction((0, 3), (2, 5)),
                         (1, 3))

    def test_tensor_contractions(self):

        with self.assertRaises(TypeError):
            generate_metric.tensor_contractions(self.M, 0)

        M_sq = sympy.Array(self.M * self.M)

        M_tensor_M = sympy.tensor.array.tensorproduct(self.M, self.M)

        # M tensor M with contraction (1,2) is exactly like a matrix product
        self.assertEqual(generate_metric.tensor_contractions(M_tensor_M,
                                                             (1, 2)), M_sq)

        expected_out = self.t**2 * self.a**2 + 2*self.t**2 + self.x**2 + self.y**2 + self.z**4

        self.assertEqual(
            generate_metric.tensor_contractions(M_tensor_M, [(0, 2), (1, 3)]),
            expected_out)

    def test_compute_Upsilon(self):
        u = self.uUP

        expected_out = -sympy.Matrix([[-self.a*u[0] + u[3], 0, 0, u[0]],
                                     [0, u[1], 0, 0],
                                     [0, 0, u[2], 0],
                                     [0, 0, 0, 2*self.z*u[3]]])

        self.assertEqual(generate_metric.compute_Upsilon(self.g_deriv, self.uUP),
                         expected_out)

    # TODO: Write test for compute_Gamma
    # def test_compute_Gamma(self):
    #     t, x, y, z, a = sympy.symbols('t, x, y, z, a', real=True)

    #     g = sympy.Matrix([[-t * a, t, 0, 0], [0, x, 0, 0], [0, 0, y, 0],
    #                       [0, 0, 0, z*z]])

    #     g_deriv = sympy.derive_by_array(g, [t, x, y, z])

    #     self.assertEqual(generate_metric.compute_Upsilon(g_deriv, uUP),
    #                      expected_out)

    def test_compute_Phi(self):
        expected_out = [sympy.Matrix([[-self.a, 0, 0, 1],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0]]),
                        sympy.Matrix([[0, 0, 0, 0],
                                      [0, 1, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0]]),
                        sympy.Matrix([[0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 0]]),
                        sympy.Matrix([[1, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 2*self.z]])]

        self.assertEqual(generate_metric.compute_Phi(self.g_deriv),
                         expected_out)

    def test_restructure_derivatives_to_list(self):

        expected_out = [sympy.Matrix([[-self.a, 0, 0, 1],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [1, 0, 0, 0]]),
                        sympy.Matrix([[0, 0, 0, 0],
                                      [0, 1, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0]]),
                        sympy.Matrix([[0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 0]]),
                        sympy.Matrix([[0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 0],
                                      [0, 0, 0, 2*self.z]])]

        self.assertEqual(generate_metric.restructure_derivatives_to_list(self.g_deriv),
                         expected_out)

    def test_generate_geodesic_clcode(self):

        with self.assertRaises(ValueError):
            generate_metric.generate_geodesic_clcode(self.g_ks_deriv,
                                                     self.g_ksUP, self.uUP, 0)

        # TODO: Improve testing. At the moment, it relies on matching the output
        #       string with an hardcoded one. This is not reliable, it is not
        #       easy to maintain, the output of the test is not informative, and
        #       tests take a long time to run

        tot_number_of_levels = 9
        level_out = {}
        for level in range(1, tot_number_of_levels + 1):
            level_out[level] = read_text('gray.etc', f'ks_level{level}.cl')
            with self.subTest(level=level):
                self.assertMultiLineEqual(generate_metric.generate_geodesic_clcode(self.g_ks_deriv,
                                                                                   self.g_ksUP,
                                                                                   self.uUP, level),
                                          level_out[level])
