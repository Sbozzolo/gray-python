#!/usr/bin/env python3

import unittest
import unittest.mock
import sympy
from sympy.codegen.rewriting import create_expand_pow_optimization
from io import StringIO
from gray import codegen


class TestCodegen(unittest.TestCase):
    def test_replace_letters_to_numbers(self):
        self.assertEqual(codegen.replace_letters_to_numbers('0123456789'),
                         'abcdefghij')
        self.assertEqual(codegen.replace_letters_to_numbers('0abc'),
                         'aabc')
        with self.assertRaises(TypeError):
            codegen.replace_letters_to_numbers(42)

    def test_lettered_symbols(self):
        xa = sympy.symbols('xa')
        self.assertEqual(next(codegen.lettered_symbols()), xa)

        prefixa = sympy.symbols('prefixa')
        self.assertEqual(next(codegen.lettered_symbols(prefix='prefix')),
                         prefixa)

        stringa = 'xa'
        self.assertEqual(next(codegen.lettered_symbols(cls=str)), stringa)

        # HACK: Not the most elegant way to test if the code actually generates
        #       the intended stream
        list_strings = ['xa', 'xb', 'xc']
        i = 0
        for e in codegen.lettered_symbols(cls=str):
            self.assertEqual(e, list_strings[i])
            i += 1
            if (i == 3):
                break

    def test_check_for_new_symbols_and_print_gray_code(self):
        t, x, y, z = sympy.symbols('t, x, y, z', real=True)

        def metric1():
            return (-t, 0, 0, 0, 0, x, 0, 0, 0, 0, y, 0, 0, 0, 0, z)

        def metric2():
            a = sympy.symbols('a', real=True)
            return (-t * a, 0, 0, 0, 0, x, 0, 0, 0, 0, y, 0, 0, 0, 0, z)

        self.assertFalse(
            codegen.check_for_new_symbols_and_print_gray_code(metric1()))

        with unittest.mock.patch('sys.stdout', new=StringIO()) as gray_code:
            expected_output = """\
Additional code in conf.c. Add after comment ADDITIONAL CODE
PARA(a, a) = atof(val);\n

Additional code in gray.c. Add after comment ADDITIONAL CODE
gray->a = 0;\n

Additional code in gray.h. Add after comment ADDITIONAL CODE
double a;\n
"""
            self.assertTrue(
                codegen.check_for_new_symbols_and_print_gray_code(
                    metric2()))

            self.assertEqual(gray_code.getvalue(), expected_output)

    def test_assign_tranposed_matrix(self):

        expected_out = "real16 UpsilonT = Upsilon.s048c159d26ae37bf;\n"
        self.assertEqual(codegen.assign_transposed_matrix("Upsilon"),
                         expected_out)

        expected_outT = "real16 Upsilon = UpsilonT.s048c159d26ae37bf;\n"
        self.assertEqual(
            codegen.assign_transposed_matrix("UpsilonT",
                                                  already_transposed=True),
            expected_outT)

        with self.assertRaises(ValueError):
            codegen.assign_transposed_matrix("Upsilon",
                                                  already_transposed=True)

    def test_find_all_numbers_in_string(self):
        with self.assertRaises(TypeError):
            codegen.find_all_numbers_in_string(0)

        self.assertEqual(codegen.find_all_numbers_in_string("123"),
                         {"123"})

        self.assertEqual(codegen.find_all_numbers_in_string("a123"),
                         {"123"})

        self.assertEqual(codegen.find_all_numbers_in_string("a1.23"),
                         {"1.23"})

        self.assertEqual(codegen.find_all_numbers_in_string("a1.2b3"),
                         {"1.2", "3"})

        self.assertEqual(codegen.find_all_numbers_in_string("a1.2b3e8"),
                         {"1.2", "3e8"})

        self.assertEqual(codegen.find_all_numbers_in_string("a1.2b3-8"),
                         {"1.2", "3", "-8"})

        self.assertEqual(codegen.find_all_numbers_in_string("a.4c.1"),
                         {".4", ".1"})

        self.assertEqual(codegen.find_all_numbers_in_string("12.e"),
                         {"12."})

        self.assertEqual(codegen.find_all_numbers_in_string("12.e4asd"),
                         {"12.e4"})

    def test_preface_numbers_with_K(self):
        with self.assertRaises(TypeError):
            codegen.preface_numbers_with_K(0)

        self.assertEqual(codegen.preface_numbers_with_K("abc"), "abc")

        self.assertEqual(codegen.preface_numbers_with_K("abc1"),
                         "abcK(1)")

        self.assertEqual(codegen.preface_numbers_with_K("1.0abc1"),
                         "K(1.0)abcK(1)")

        self.assertEqual(codegen.preface_numbers_with_K("1.0abc2.1"),
                         "K(1.0)abcK(2.1)")

        self.assertEqual(codegen.preface_numbers_with_K("1abc2.1"),
                         "K(1)abcK(2.1)")

        self.assertEqual(codegen.preface_numbers_with_K("1abc2.1efg2.2"),
                         "K(1)abcK(2.1)efgK(2.2)")

        self.assertEqual(
            codegen.preface_numbers_with_K("1abc2.1efg2.2i1e2"),
            "K(1)abcK(2.1)efgK(2.2)iK(1e2)")

        self.assertEqual(codegen.preface_numbers_with_K("1abc/2.1"),
                         "K(1)abc/K(2.1)")

    def test_replace_powers_with_gray_functions(self):

        gray_square = sympy.Function('GRAY_SQUARE')
        gray_cube = sympy.Function('GRAY_CUBE')
        gray_four = sympy.Function('GRAY_FOUR')
        gray_sqrt = sympy.Function('GRAY_SQRT')
        gray_sqrt_cube = sympy.Function('GRAY_SQRT_CUBE')

        repl = {
            2: gray_square,
            3: gray_cube,
            4: gray_four,
            0.5: gray_sqrt,
            1.5: gray_sqrt_cube
        }

        a, b, c = sympy.symbols('a b c')

        expand_opt4 = create_expand_pow_optimization(4)

        with self.assertRaises(TypeError):
            codegen.replace_powers_with_gray_functions(0, repl)

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(a**2, repl),
            expand_opt4(a * a))

        self.assertEqual(
            codegen.replace_powers_with_gray_functions((a + b)**2, repl),
            gray_square(a + b))

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(
                1 / (a + b)**2, repl), 1 / gray_square(a + b))

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(
                1 / (a + b)**2 + a**5, repl), 1 / gray_square(a + b) + a**5)

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(
                1 / (a + b)**2 + a**5 + sympy.sqrt(b), repl),
            1 / gray_square(a + b) + a**5 + gray_sqrt(b))

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(
                1 / (a + b)**2 + a**5 + sympy.sqrt(a + b + c), repl),
            1 / gray_square(a + b) + a**5 + gray_sqrt(a + b + c))

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(
                1 / (a + b)**(3 / 2), repl), 1 / gray_sqrt_cube(a + b))

        self.assertEqual(
            codegen.replace_powers_with_gray_functions(
                1 / (a + b)**(1.5), repl), 1 / gray_sqrt_cube(a + b))

    def test_convert_expr_to_C(self):
        with self.assertRaises(TypeError):
            codegen.convert_expr_to_C(0)

        a, b, c = sympy.symbols('a b c')

        # ccode changes the order of factors, this is why some asserts look
        # weird. What is important is that it is mathematically correct

        self.assertEqual(codegen.convert_expr_to_C(sympy.Float(2)),
                         "K(2.0)")

        self.assertEqual(codegen.convert_expr_to_C(2.0 * a), "K(2.0)*a")

        self.assertEqual(codegen.convert_expr_to_C(2.0 * a * a),
                         "K(2.0)*(a*a)")

        self.assertEqual(codegen.convert_expr_to_C(2.0 * a**5.0),
                         "K(2.0)*pow(a, K(5.0))")

        self.assertEqual(codegen.convert_expr_to_C(sympy.sqrt(b) * a),
                         "a*GRAY_SQRT(b)")

        self.assertEqual(codegen.convert_expr_to_C(1.0 / (a + b)**2),
                         "K(1.0)/GRAY_SQUARE(a + b)")

        self.assertEqual(codegen.convert_expr_to_C(a + b), "a + b")

        self.assertEqual(codegen.convert_expr_to_C(a + b + 1.0),
                         "a + b + K(1.0)")

        self.assertEqual(codegen.convert_expr_to_C(a + b + 1.0),
                         "a + b + K(1.0)")

        self.assertEqual(codegen.convert_expr_to_C((a + b + c)**(3 / 2)),
                         "GRAY_SQRT_CUBE(a + b + c)")

    def test_assign_scalar(self):

        a, b = sympy.symbols('a b')

        with self.assertRaises(TypeError):
            codegen.assign_scalar(0, a)

        with self.assertRaises(TypeError):
            codegen.assign_scalar('a', 0)

        self.assertEqual(codegen.assign_scalar('a_name', a * a),
                         "real a_name = a*a;\n")

        self.assertEqual(codegen.assign_scalar('a_name', (a + b)**2),
                         "real a_name = GRAY_SQUARE(a + b);\n")

    def test_assign_vector(self):
        a, b, c, d = sympy.symbols('a b c d')

        with self.assertRaises(TypeError):
            codegen.assign_vector(0, [a])

        with self.assertRaises(TypeError):
            codegen.assign_vector('a', 0)

        self.assertEqual(codegen.assign_vector('v_name', [a, b, c, d]),
                         "real4 v_name = {a,\nb,\nc,\nd};\n")

    def test_assign_vectors(self):
        a, b, c, d = sympy.symbols('a b c d')

        with self.assertRaises(TypeError):
            codegen.assign_vectors(0)

        expected_out = "real4 v_name = {a,\nb,\nc,\nd};\n"
        expected_out += "real4 v_name2 = {a,\nb,\nd,\nc};\n"

        vector_dict = {'v_name': [a, b, c, d], 'v_name2': [a, b, d, c]}

        self.assertEqual(codegen.assign_vectors(vector_dict),
                         expected_out)

    def test_assign_matrix(self):
        a, b, c, d = sympy.symbols('a b c d')
        e, f, g, h = sympy.symbols('e f g h')
        i, l, m, n = sympy.symbols('i l m n')
        o, p, q, r = sympy.symbols('o p q r')

        M = sympy.Matrix([[a, b, c, d], [e, f, g, h], [i, l, m, n],
                          [o, p, q, r]])

        expected_out = "real16 m_name = \
{a,\nb,\nc,\nd,\ne,\nf,\ng,\nh,\ni,\nl,\nm,\nn,\no,\np,\nq,\nr};\n"

        with self.assertRaises(TypeError):
            codegen.assign_matrix(0, M)

        with self.assertRaises(TypeError):
            codegen.assign_matrix('a', 0)

        self.assertEqual(codegen.assign_matrix('m_name', M), expected_out)

    def test_assign_matrices(self):
        a, b, c, d = sympy.symbols('a b c d')
        e, f, g, h = sympy.symbols('e f g h')
        i, l, m, n = sympy.symbols('i l m n')
        o, p, q, r = sympy.symbols('o p q r')

        M = sympy.Matrix([[a, b, c, d], [e, f, g, h], [i, l, m, n],
                          [o, p, q, r]])

        expected_out = "real16 m_name = \
{a,\nb,\nc,\nd,\ne,\nf,\ng,\nh,\ni,\nl,\nm,\nn,\no,\np,\nq,\nr};\n"
        expected_out += "real16 m_name2 = \
{a,\nb,\nc,\nd,\ne,\nf,\ng,\nh,\ni,\nl,\nm,\nn,\no,\np,\nq,\nr};\n"

        with self.assertRaises(TypeError):
            codegen.assign_matrices(0)

        M_dict = {'m_name': M, 'm_name2': M}

        self.assertEqual(codegen.assign_matrices(M_dict), expected_out)

    def test_apply_cse_and_assign_symbols(self):

        a, b = sympy.symbols('a b')
        gray_suba = sympy.symbols('gray_suba')

        expected_out = [gray_suba**2 + gray_suba], "real gray_suba = a + b;\n"

        self.assertEqual(
            codegen.apply_cse_and_assign_symbols(a + b + (a + b)**2),
            expected_out)

        self.assertEqual(
            codegen.apply_cse_and_assign_symbols([a + b + (a + b)**2]),
            expected_out)

        expected_out = [gray_suba**2 + gray_suba,
                        gray_suba], "real gray_suba = a + b;\n"

        self.assertEqual(
            codegen.apply_cse_and_assign_symbols(
                [a + b + (a + b)**2, a + b]), expected_out)
