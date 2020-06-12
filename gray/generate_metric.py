#!/usr/bin/env python3

"""gray_generate_metric.py provides functions to produce OpenCL kernels for arbitrary
spacetimes.

"""

from sympy import tensorproduct, tensorcontraction
from sympy import Matrix, Array
from gray.codegen import apply_cse_and_assign_symbols
from gray.codegen import assign_transposed_matrix
from gray.codegen import assign_scalar, assign_vector, assign_matrix
from gray.codegen import assign_matrices, assign_vectors


def metric_from_tuple(tup):
    """Return a sympy Matrix from a 16-elements tuple"""

    if (not type(tup) is tuple):
        raise TypeError("tup in metric_from_tuple has to be a tuple")

    if (not len(tup) == 16):
        raise ValueError("tup metric_from_tuple has to have 16 elements")

    (g_tt, g_tx, g_ty, g_tz, _, g_xx, g_xy, g_xz, _, _, g_yy, g_yz, _, _, _,
     g_zz) = tup

    return Matrix([[g_tt, g_tx, g_ty, g_tz], [g_tx, g_xx, g_xy, g_xz],
                   [g_ty, g_xy, g_yy, g_yz], [g_tz, g_xz, g_yz, g_zz]])


def find_indeces_for_contraction(previous_indeces, current_indeces):
    """Sympy only supports contraction of one index at the time.

To have multiple contractions, we have to nest tensorcontraction calls. The
problem is that each time the indeces structure changes. Consider for example
the contraction

A^abc B_cda (0)

Python indeces are

A^012 B_345

So, we would like to say "contract [(0,5), (2,3)]" to achieve (0).
This is done by first contracting the first, which results in

A^bc B_cd

which has indeces

A^01 B_23

So the next contraction has to be (1,2), not (2, 3).

This function takes care of adjusting the indeces.

    """
    # Tuples do not support assignment,
    # we revert back to tuple at the return
    new_indeces = list(current_indeces)

    for i, cur_ind in enumerate(current_indeces):
        for prev_ind in previous_indeces:
            if cur_ind > prev_ind:
                new_indeces[i] -= 1
    return tuple(new_indeces)


def tensor_contractions(T, indeces):
    """If indeces is a tuple with two integers, then perform the contraction
    between those two indeces. If indeces is a list of tuples, then apply
    multiple contractions.
    """
    if (type(indeces) is tuple):
        ret = tensorcontraction(T, indeces)
    elif (type(indeces) is list):
        T = tensorcontraction(T, indeces[0])
        for i in range(1, len(indeces)):
            T = tensorcontraction(T, find_indeces_for_contraction(indeces[i-1],
                                                                  indeces[i]))
        ret = T
    else:
        raise TypeError("indeces in tensor_contractions has to be a list of \
tuples or a tuple")

    return ret


def compute_Upsilon(g_deriv, uUP):
    """Return Upsilon =

    Upsilon_ alpha beta = g_beta gamma, alpha uUP^gamma

    """

    gdU = tensorproduct(g_deriv, uUP)

    return -Matrix(tensorcontraction(gdU, (1, 3)))


def compute_Gamma(g_deriv, gUP):
    """Return Christoffel symbols
    """
    g_derivT = Array([(g_deriv[:, :, i]).transpose() for i in range(4)])

    gUgd = tensorproduct(gUP, g_deriv)
    gUgdT = tensorproduct(gUP, g_derivT)

    return 1 / 2 * (tensorcontraction(gUgd, (1, 3)) +
                    tensorcontraction(gUgdT, (1, 3)) -
                    tensorcontraction(gUgd, (1, 2)))


def restructure_derivatives_to_list(g_deriv):
    """Sympy returns derivative of the metric in an unfriendly format
that is not easy to use. Here, we transform it into a list of matrices,
each entry of the list is derivative with respect to a different variables.
So, for example, return[0] is derivative with respect of t.
    """

    g_deriv_t = (Matrix(g_deriv[0, :, :]))
    g_deriv_x = (Matrix(g_deriv[1, :, :]))
    g_deriv_y = (Matrix(g_deriv[2, :, :]))
    g_deriv_z = (Matrix(g_deriv[3, :, :]))

    return [g_deriv_t, g_deriv_x, g_deriv_y, g_deriv_z]


def compute_Phi(g_deriv):
    """
    Phi^(mu)_alpha beta = g_beta (mu), alpha
    """

    # Remember, the first index is the derivative index
    Phi_t = Matrix(g_deriv[:, 0, :])
    Phi_x = Matrix(g_deriv[:, 1, :])
    Phi_y = Matrix(g_deriv[:, 2, :])
    Phi_z = Matrix(g_deriv[:, 3, :])

    return [Phi_t, Phi_x, Phi_y, Phi_z]


def generate_geodesic_clcode(g_deriv, gUP, uUP, level):
    """Reduce RHS to 1 equation"""

    # TODO: This function is a little bit too complicated (cyclomatic complexity
    #       is 19...). It may be useful to break it down into smaller pieces, so
    #       that testing can also be improved.

    # in our template we have to insert:
    # - The CSE symbols
    # - Additional definitions/manipulations
    # - The RHS

    # To do this we initialize the variable clcode with the cse symbols, and we
    # append everything else to it

    if (level == 1):
        # Upsilon_ beta alpha = g_beta gamma, alpha uUP gamma
        Upsilon = compute_Upsilon(g_deriv, uUP)

        gUPuUPUpsilon = tensorproduct(tensorproduct(gUP, uUP), Upsilon)

        # Xi1^mu = g^mu alpha uUP beta Upsilon_beta alpha
        # Xi2^mu = - 1/2 g^mu alpha u^beta Upsilon_alpha beta
        Xi1UP = tensor_contractions(gUPuUPUpsilon, [(2, 3), (1, 4)])
        Xi2UP = tensor_contractions(-1 / 2 * gUPuUPUpsilon, [(1, 3), (2, 3)])

        # rhs^mu = -(Xi1^ mu + Xi2^ mu)
        rhsUP, clcode = apply_cse_and_assign_symbols(list(Xi1UP + Xi2UP))
        clcode += assign_vector('rhs', rhsUP)
    elif (level == 2):
        Upsilon = compute_Upsilon(g_deriv, uUP)
        Xi = Upsilon.transpose() - Upsilon / 2

        # Lambda = gUP * Xi
        Lambda, clcode = apply_cse_and_assign_symbols(gUP * Xi)
        clcode += assign_matrix('Lambda', Lambda[0])
        clcode += "real4 rhs = matrix_vector_product(Lambda, u);\n"
    elif (level == 3):
        Upsilon = compute_Upsilon(g_deriv, uUP)
        Xi = Upsilon.transpose() - Upsilon / 2
        (Xi, gUP), clcode = apply_cse_and_assign_symbols([Xi, gUP])

        clcode += assign_matrices({'Xi': Xi, 'gUP': gUP})
        clcode += "real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));\n"
    elif (level == 4):
        Upsilon = compute_Upsilon(g_deriv, uUP)
        (Upsilon, gUP), clcode = apply_cse_and_assign_symbols([Upsilon, gUP])

        clcode += assign_matrices({'Upsilon': Upsilon, 'gUP': gUP})
        clcode += assign_transposed_matrix('Upsilon')
        clcode += "real16 Xi = UpsilonT - Upsilon/2;\n"
        clcode += "real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));\n"
    elif (level == 5):
        Phi = compute_Phi(g_deriv)
        (*Phi, gUP), clcode = apply_cse_and_assign_symbols([*Phi, gUP])
        clcode += assign_matrices({'Phi_t': Phi[0],
                                   'Phi_x': Phi[1],
                                   'Phi_y': Phi[2],
                                   'Phi_z': Phi[3],
                                   'gUP': gUP})

        clcode += "real16 Upsilon = -uUPt * Phi_t - uUPx * Phi_x"
        clcode += " - uUPy * Phi_y - uUPz * Phi_z;\n"
        clcode += assign_transposed_matrix('Upsilon')
        clcode += "real16 Xi = UpsilonT - Upsilon/2;\n"
        clcode += "real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));\n"
    elif (level == 6):
        # Typically this level is very slow because it has to parse a lot of
        # long sympy expressions

        # CSE does not handle very well g_deriv, so we split it
        g_deriv = restructure_derivatives_to_list(g_deriv)

        (*g_deriv, gUP), clcode = apply_cse_and_assign_symbols([*g_deriv, gUP])

        d = {0: 't', 1: 'x', 2: 'y', 3: 'z'}

        # Write derivatives
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    clcode += assign_scalar("g_" + d[a] + d[b] + d[c],
                                            g_deriv[c][a, b])

        for c in range(4):
            clcode += "real16 Phi_" + d[c] + " = {" + \
                ", ".join(["g_" + d[b] + d[c] + d[a] for a in range(4)
                           for b in range(4)]) + "};\n"

        clcode += "real16 Upsilon = -uUPt * Phi_t - uUPx * Phi_x"
        clcode += " - uUPy * Phi_y - uUPz * Phi_z;\n"
        clcode += assign_transposed_matrix('Upsilon')
        clcode += assign_matrix('gUP', gUP)
        clcode += "real16 Xi = UpsilonT - Upsilon/2;\n"
        clcode += "real4 rhs = matrix_vector_product(gUP, matrix_vector_product(Xi, u));\n"
    elif (level == 7):
        Gamma = compute_Gamma(g_deriv, gUP)
        GammauUPuUP = tensorproduct(tensorproduct(Gamma, uUP), uUP)
        rhsUP = -tensor_contractions(GammauUPuUP, [(1, 3), (2, 4)])
        rhsUP, clcode = apply_cse_and_assign_symbols(list(rhsUP))
        clcode += assign_vector('rhs', rhsUP)
    elif (level == 8):
        Gamma = [Matrix(G) for G in compute_Gamma(g_deriv, gUP)]
        GammaUP, clcode = apply_cse_and_assign_symbols(Gamma)

        d = {0: 't', 1: 'x', 2: 'y', 3: 'z'}

        clcode += assign_matrices({"GammaUP" + d[i]: GammaUP[i] for i in range(4)})

        clcode += "real4 rhs = {-dot(u, matrix_vector_product(GammaUPt, u)),"
        clcode += "-dot(u, matrix_vector_product(GammaUPx, u)),"
        clcode += "-dot(u, matrix_vector_product(GammaUPy, u)),"
        clcode += "-dot(u, matrix_vector_product(GammaUPz, u))};\n"
    elif (level == 9):
        # CSE does not handle very well g_deriv, so we split it
        # and we make a list with all the derivatives
        # g_deriv[0] is derivative with respect to t, and so on
        g_deriv = restructure_derivatives_to_list(g_deriv)

        # We unpack to avoid nested lists
        (*g_deriv, gUP), clcode = apply_cse_and_assign_symbols([*g_deriv, gUP])

        d = {0: 't', 1: 'x', 2: 'y', 3: 'z'}

        # Write derivatives
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    clcode += assign_scalar("g_" + d[a] + d[b] + d[c],
                                            g_deriv[c][a, b])

        clcode += assign_vectors({'gUP' + d[i]: list(gUP.row(i))
                                  for i in range(4)})

        for c in range(4):
            clcode += "real16 A_" + d[c] + " = {" + \
                ", ".join(["g_" + d[c] + d[a] + d[b] for a in range(4)
                           for b in range(4)]) + "};\n"
            clcode += "real16 B_" + d[c] + " = {" + \
                ", ".join(["g_" + d[c] + d[b] + d[a] for a in range(4)
                           for b in range(4)]) + "};\n"
            clcode += "real16 C_" + d[c] + " = {" + \
                ", ".join(["g_" + d[a] + d[b] + d[c] for a in range(4)
                           for b in range(4)]) + "};\n"
        for c in range(4):
            clcode += "real16 GammaUP" + d[c] + " = 0.5 * (" + \
                "gUP" + d[c] + ".s0 * (A_t + B_t - C_t) + " \
                "gUP" + d[c] + ".s1 * (A_x + B_x - C_x) + " \
                "gUP" + d[c] + ".s2 * (A_y + B_y - C_y) + " \
                "gUP" + d[c] + ".s3 * (A_z + B_z - C_z));\n"

        clcode += "real4 rhs = {-dot(u, matrix_vector_product(GammaUPt, u)),"
        clcode += "-dot(u, matrix_vector_product(GammaUPx, u)),"
        clcode += "-dot(u, matrix_vector_product(GammaUPy, u)),"
        clcode += "-dot(u, matrix_vector_product(GammaUPz, u))};\n"
    else:
        raise ValueError("Level {} not implemented".format(level))

    return clcode
