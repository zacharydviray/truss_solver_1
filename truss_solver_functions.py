import numpy as np

import math


def assign_BCs(NL, ENL):
    PD = np.size(NL, 1)
    NoN = np.size(NL, 0)

    DOFs = 0
    DOCs = 0

    # Assign DOF / constrained DOF numbering
    for i in range(NoN):
        for j in range(PD):

            if ENL[i, PD + j] == -1:
                DOCs -= 1
                ENL[i, 2 * PD + j] = DOCs
            else:
                DOFs += 1
                ENL[i, 2 * PD + j] = DOFs

    # Map constrained DOFs to global numbering
    for i in range(NoN):
        for j in range(PD):

            if ENL[i, 2 * PD + j] < 0:
                ENL[i, 3 * PD + j] = abs(ENL[i, 2 * PD + j]) + DOFs
            else:
                ENL[i, 3 * PD + j] = ENL[i, 2 * PD + j]

    DOCs = abs(DOCs)

    return ENL, DOFs, DOCs


# ----------------------------
# GLOBAL STIFFNESS ASSEMBLY
# ----------------------------
def assemble_stiffness(ENL, EL, NL, E, A):

    NoE = np.size(EL, 0)
    NPE = np.size(EL, 1)
    PD = np.size(NL, 1)
    NoN = np.size(NL, 0)

    K = np.zeros((NoN * PD, NoN * PD))

    for i in range(NoE):

        n1 = EL[i, 0:NPE].astype(int)

        k = element_stiffness(n1, ENL, E, A)

        for r in range(NPE):
            for p in range(PD):
                for q in range(NPE):
                    for s in range(PD):

                        row = int(ENL[n1[r] - 1, 3 * PD + p]) - 1
                        col = int(ENL[n1[q] - 1, 3 * PD + s]) - 1

                        K[row, col] += k[r * PD + p, q * PD + s]

    return K


# ----------------------------
# ELEMENT STIFFNESS MATRIX
# ----------------------------
def element_stiffness(n1, ENL, E, A):

    X1 = ENL[n1[0] - 1, 0]
    Y1 = ENL[n1[0] - 1, 1]

    X2 = ENL[n1[1] - 1, 0]
    Y2 = ENL[n1[1] - 1, 1]

    L = math.sqrt((X2 - X1) ** 2 + (Y2 - Y1) ** 2)

    if np.isclose(L, 0):
        raise ValueError(f"Zero-length element detected: {n1}")

    C = (X2 - X1) / L
    S = (Y2 - Y1) / L

    k = (E * A / L) * np.array([
        [C * C, C * S, -C * C, -C * S],
        [C * S, S * S, -C * S, -S * S],
        [-C * C, -C * S, C * C, C * S],
        [-C * S, -S * S, C * S, S * S]
    ])

    return k


def assemble_forces(ENL, NL):
    PD = np.size(NL, 1)
    NoN = np.size(NL, 0)

    Fp = []

    for i in range(NoN):
        for j in range(PD):
            if ENL[i, PD + j] == 1:
                Fp.append(ENL[i, 5 * PD + j])

    Fp = np.vstack(Fp).reshape(-1, 1)

    return Fp


def assemble_displacement(ENL, NL):
    PD = np.size(NL, 1)
    NoN = np.size(NL, 0)

    Up = []

    for i in range(NoN):
        for j in range(PD):
            if ENL[i, PD + j] == -1:
                Up.append(ENL[i, 5 * PD + j])

    Up = np.vstack(Up).reshape(-1, 1)

    return Up


def update_nodes(ENL, U_u, NL, Fu):
    PD = np.size(NL, 1)
    NoN = np.size(NL, 0)

    free_counter = 0
    constrained_counter = 0

    for i in range(NoN):
        for j in range(PD):

            if ENL[i, PD + j] == 1:  # free DOF
                free_counter += 1
                ENL[i, 4 * PD + j] = U_u[free_counter - 1, 0]

            else:  # constrained DOF
                constrained_counter += 1
                ENL[i, 5 * PD + j] = Fu[constrained_counter - 1, 0]

    return ENL























