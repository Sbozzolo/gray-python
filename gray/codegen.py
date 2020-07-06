"""gray_codegen.py provides functions to produce OpenCL kernels.

"""

import re
from sympy import symbols, Matrix, Basic, Wild, Function, cse
from sympy import ImmutableDenseMatrix, MutableDenseMatrix
from sympy.codegen.rewriting import create_expand_pow_optimization
from string import ascii_letters
from sympy.printing import ccode


def replace_letters_to_numbers(text):
    """Replace all the integers in text with the letters in the alphabet,
    for example: 0 -> a, 1 -> b, and so on.
    """

    if (type(text) != str):
        raise TypeError("replace_letters_to_numbers can only handle strings!")

    # We have to do some manipulations with numbers,
    # so it is best if the symbol names have no numbers
    # but letters, so that they don't clash
    # This is essentially
    # '0': 'a',
    # '1': 'b',
    # .......
    # '9': 'j'
    replacements = {str(i): ascii_letters[i] for i in range(10)}

    return re.sub(r'\d', lambda m: replacements[m.group()], text)


def lettered_symbols(prefix='x', cls=None, start=0):
    """Generate an infinite stream of Symbols consisting of a prefix and
    increasing alphabetical subscripts.

    """
    if cls is None:
        # We can't just make the default cls=Symbol because it isn't
        # imported yet.
        from sympy import Symbol
        cls = Symbol
    while True:
        name = '%s%s' % (prefix, replace_letters_to_numbers(str(start)))
        s = cls(name)
        yield s
        start += 1


def assign_transposed_matrix(matrix_name, already_transposed=False):
    """Print OpenCL code to define a matrix matrix_nameT that is
    the transpose of the matrix.

    If already_transposed is True, then defined matrix_name (without the T)

    """

    if (already_transposed):
        if (not matrix_name.endswith("T")):
            raise ValueError("Matrix name does not suggest that the matrix is\
                              already transposed")
        # Remove "T"
        transposed_name = matrix_name[:-1]
    else:
        transposed_name = matrix_name + "T"

    # The way we access matrices in OpenCL is with the indeces
    # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f
    # This is the matrix
    #    0 1 2 3
    #    4 5 6 6
    #    8 9 a b
    #    c d e f
    # So the transposed is
    #    0 4 8 c
    #    1 5 9 d
    #    2 6 a e
    #    3 6 b f

    return f"real16 {transposed_name} = {matrix_name}.s048c159d26ae37bf;\n"


def check_for_new_symbols_and_print_gray_code(metric):
    """Check if metric() has any other symbol except t, x, y, z. In that case,
    print the code that has to be added to GRay2 to use those symbols.

    The argument metric is a tuple with the 16 components of the metric.

    Additional parameters are always initialized with the default value of 0.
    If you want to change that, change the lines that go in the file gray.c

    If new symbols are found, return True, otherwise return False.

    """

    new_symbols = False
    known_symbols = set(symbols('t, x, y, z', real=True))

    # File param.opts
    param_code = ""

    # Note that .free_symbols only returns free symbols. For example, for
    # Sum(T, (n, 1, N))/N it returns {N, T}, but not n. This is what we want.
    # Note also that the output of .free_symbols is a Python set.

    # The output of metric is a tuple, so we convert it to a sympy Matrix
    symbols_in_metric = Matrix(metric).free_symbols

    # Unknown symbols are the difference between the two sets
    unkown_symbols = symbols_in_metric - known_symbols

    # Now we print the code
    for s in unkown_symbols:
        new_symbols = True
        str_s = str(s)
        param_code += f"double {str_s}:{str_s} = strtod(val, &rem);\n"

    if (new_symbols):
        print("Additional code in param.opts:")
        print(param_code)

    return (new_symbols)


def find_all_numbers_in_string(text):
    """Find all the numbers in the string text and return a set with all
    the unique occurrences.

    """
    if (type(text) != str):
        raise TypeError("find_all_numbers_in_string can only handle strings!")

    # From: StackOverflow
    # (4703390/how-to-extract-a-floating-number-from-a-string)
    numeric_const_pattern = r"""
        [-+]? # optional sign
        (?:
            (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
            |
            (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
        )
        (?: [Ee] [+-]? \d+ ) ?
    """

    # Compile the pattern
    rx = re.compile(numeric_const_pattern, re.VERBOSE)

    # Here we use set() for to get the unique elements
    return set(rx.findall(text))


def preface_numbers_with_K(text):
    """To make sure that floats have consistent precision, in GRay, we pass floats
to the macro K(x). This function finds all the numbers in text and substitute
them whit K(number).

    """
    if (type(text) != str):
        raise TypeError("preface_numbers_with_K can only handle strings!")

    subs = find_all_numbers_in_string(text)
    if (len(subs) > 0):
        # find_all_numbers_in_string returns a set. We transform back to a list
        # so that we can order it
        subs = list(subs)
        # The reason is that we may have to do substitutions like
        # 1.0 and 1.00023
        # We first have to substitute the second, otherwise it will
        # become K(1.0)00023
        subs = sorted(subs, key=len)
        # Long elements first
        subs = reversed(subs)

        subs_dict = {s: f'K({s})' for s in subs}

        # Substitute everything
        # From StackOverflow
        # (15175142/how-can-i-do-multiple-substitutions-using-regex-in-python)

        # Create a regular expression from the dictionary keys
        regex = re.compile("(%s)" % "|".join(map(re.escape, subs_dict.keys())))

        # For each match, look-up corresponding value in dictionary
        ret = regex.sub(lambda mo: subs_dict[mo.string[mo.start():mo.end()]],
                         text)
    else:
        ret = text

    return ret


def replace_powers_with_gray_functions(expr, repl, up_to_exponent=4):
    """repls is a dictionary of the form:
    {exponent_of_the_power: sympy.Function("GRAY_FUNCTION")}.

    GRAY_FUNCTIONS have to be defined in GRay! The definition can be added
    to the template of the output file.

    """
    if (not isinstance(expr, Basic)):
        raise TypeError(
            "replace_powers_with_gray_functions can only handle sympy expressions!"
        )

    # Instead of having pow(x,3), write x*x*x
    # create_expand_pow_optimization(N) performs this substitution
    # up to the power N
    expand_opt = create_expand_pow_optimization(up_to_exponent)

    expr = expand_opt(expr)

    # The previous substitution does not work when there are multiple
    # symbols involved, for example (a+b)**2
    # Having performed cse, these expressions should not be too bad

    a, b = Wild('a'), Wild('b')

    for e in repl.keys():
        expr = expr.replace((a + b)**e, lambda a, b: repl[e](a + b))
        expr = expr.replace((a + b)**-e, lambda a, b: 1 / repl[e](a + b))

    # For the square roots we also match single symbols
    # Square roots are identified by non-integer exponents
    sqrt_elems = [e for e in repl.keys() if type(e) is float]
    for e in sqrt_elems:
        expr = expr.replace(a**e, lambda a: repl[e](a))
        expr = expr.replace(a**-e, lambda a: 1 / repl[e](a))

    return expr


def convert_expr_to_C(expr):
    """Convert expression expr to C (OpenCL) code and convert powers to
    GRay functions, like (a+b)**2 to GRAY_SQUARE(a+b). The conversions
    currently performed are:

        x**2: GRAY_SQUARE(x)
        x**3: GRAY_CUBE(x)
        x**4: GRAY_FOUR(x)
        x**0.5: GRAY_SQRT(x)
        x**1.5: GRAY_SQRT_CUBE(x)

    GRAY_FUNCTIONS are defined in GRay through the file template.
    """

    if (not isinstance(expr, Basic)):
        raise TypeError("convert_expr_to_C can only handle sympy expressions!")

    # To work around this, we pattern-match all those expressions and
    # manually substitute them with a temporary function GRAY_SQUARE, which
    # is defined in GRay and takes the square. We do the same for cubes,
    # square roots, forth powers, and square roots of cubes.
    # HACK: It would be best to find a way to do this directly here, and
    #       hot in GRay

    # NOTE: These functions (GRAY_SQUARE, etc) are defined in GRay, if you
    #       add additional functions, make sure to also add the definition
    repl = {
        2: Function('GRAY_SQUARE'),
        3: Function('GRAY_CUBE'),
        4: Function('GRAY_FOUR'),
        0.5: Function('GRAY_SQRT'),
        1.5: Function('GRAY_SQRT_CUBE')
    }

    expr = replace_powers_with_gray_functions(expr, repl)

    # The sympy expression contains the functions GRAY_SQUARE and friends, but
    # sympy does not know how to translate them to C code. We tell sympy to
    # translate GRAY_SQUARE with GRAY_SQUARE (and so on) by providing the
    # user_functions dictionary
    user_functions = {str(v): str(v) for v in repl.values()}
    string_expr = ccode(expr, user_functions=user_functions)

    # GRay handles precision of floats in a very peculiar way.

    # To make sure that floats have consistent precision, in GRay, we pass
    # floats to the macro K(x).
    # Here, we have to do the same. We use a regular expression to match
    # all the numbers, then, we substitute them with K(NUMBER), where
    # NUMBER is the match float
    return preface_numbers_with_K(string_expr)


def assign_scalar(scalar_name, scalar):
    """Return OpenCL code to assign a scalar (a real)

    scalar_name is the name of variable in OpenCL, scalar the sympy expression

    """
    if (not isinstance(scalar, Basic)):
        raise TypeError("scalar in assign_scalar has to be a sympy expression")

    if (not type(scalar_name) is str):
        raise TypeError("scalar_name in assign_scalar has to be string")

    return f"real {scalar_name} = {convert_expr_to_C(scalar)};\n"


def assign_vector(vector_name, vector):
    """Return OpenCL code to assign a vector (a real4)

    vector_name is the name of variable in OpenCL, vector is a list of
    sympy expressions

    """
    if (not type(vector) is list):
        raise TypeError("vector in assign_vector has to be \
a list of sympy expressions")

    if (not type(vector_name) is str):
        raise TypeError("vector_name in assign_vector has to be string")

    assign = "real4 {0} = {{{1},\n{2},\n{3},\n{4}}};\n".format(
        vector_name, *list(map(convert_expr_to_C, vector)))
    return assign


def assign_matrix(matrix_name, matrix):
    """Return OpenCL code to assign a matrix (a real16)

    matrix_name is the name of variable in OpenCL, matrix is sympy
    Matrix or a list
    """

    if (not (isinstance(matrix, (ImmutableDenseMatrix, MutableDenseMatrix)))):
        raise TypeError("matrix in assign_matrix has to be a sympy Matrix")

    if (not type(matrix_name) is str):
        raise TypeError("matrix_name in assign_matrix has to be string")

    # TODO: This is very ugly looking...
    assign = "real16 {0} \
= {{{1},\n{2},\n{3},\n{4},\
\n{5},\n{6},\n{7},\n{8},\
\n{9},\n{10},\n{11},\n{12},\
\n{13},\n{14},\n{15},\n{16}}};\n".format(matrix_name,
                                         *list(map(convert_expr_to_C, matrix)))

    return assign


def assign_matrices(matrices_dictionary):
    """Return OpenCL code to assign multiple matrices (real16).

Input is a dictionary of the form name: sympy_expression

    """
    if (not type(matrices_dictionary) is dict):
        raise TypeError("matrix_dictionary in assign_matrices has to be dictionary")

    assign = ""
    for k, v in matrices_dictionary.items():
        assign += assign_matrix(k, v)
    return assign


def assign_vectors(vectors_dictionary):
    """Return OpenCL code to assign multiple vectors (real4).

Input is a dictionary of the form {name: list of sympy expressions}

    """
    if (not type(vectors_dictionary) is dict):
        raise TypeError("vectors_dictionary in assign_vectors has to be dictionary")

    assign = ""
    for k, v in vectors_dictionary.items():
        assign += assign_vector(k, v)
    return assign


def apply_cse_and_assign_symbols(exprs):
    """We stand at the crossroads: what is more convenient, to our expressions
    as functions of the basic symbols (t,x,y,z), possibly allowing for greatest
    symbolical simplification; or to introduce temporary symbols in case
    simplifications are not possible?

    In very simple cases, sympy can reduce expressions to their most compact
    and computationally cheap form. In this case, we can have a single set of
    equations for the right-hand-side of the geodesic equations, and this would
    be most likely the best way to write these equations. However, reality is
    often disappointing, and sympy cannot reduce effectively equations, or the
    results are terribly long and complicated, with too many operations. In
    this latter case, it is better to introduce temporary symbols for the
    common factors, performing Common Sub-expression Elimination (CSE).

    Sympy has a function, cse, to list all the recognized sub-expressions.
    This routine returns a list of pairs symbols and expressions.

    apply_cse_and_assign_sybols applies cse to the list exprs returns the
    reduced expression as a function of those symbols, and the OpenCL code.

    Regardless of the input, the reduced expression(s) are always in a list.

    """

    repls, reduced = cse(exprs, symbols=lettered_symbols(prefix='gray_sub'))
    cse_symbols = ""
    for k, v in repls:
        cse_symbols += assign_scalar(str(k), v)
    return reduced, cse_symbols
