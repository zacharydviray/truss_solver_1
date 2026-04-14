import numpy as np
from truss_solver_functions import *

# -------------------------
# INPUT
# -------------------------
NL = np.array([[0, 0],
               [1, 0],
               [0.5, 1]])

EL = np.array([[1, 2],
               [2, 3],
               [3, 1]])

DorN = np.array([[-1, -1],
                 [1, -1],
                 [1, 1]])

F_input = np.array([[0, 0],
                    [0, 0],
                    [0, -20]])

U_input = np.array([[0, 0],
                    [0, 0],
                    [0, 0]])

E = 1e6
A = 0.01

PD = NL.shape[1]
NoN = NL.shape[0]

# -------------------------
# INITIALIZE ENL
# -------------------------
ENL = np.zeros((NoN, 6 * PD))

ENL[:, 0:PD] = NL
ENL[:, PD:2*PD] = DorN

# -------------------------
# DOF ASSIGNMENT
# -------------------------
ENL, DOFs, DOCs = assign_BCs(NL, ENL)

# -------------------------
# LOAD INPUTS INTO ENL
# -------------------------
ENL[:, 4*PD:5*PD] = U_input   # prescribed displacements
ENL[:, 5*PD:6*PD] = F_input   # applied forces

# -------------------------
# STIFFNESS MATRIX
# -------------------------
K = assemble_stiffness(ENL, EL, NL, E, A)

# -------------------------
# ASSEMBLE VECTORS
# -------------------------
Fp = assemble_forces(ENL, NL)        # forces at free DOFs
Up = assemble_displacement(ENL, NL)  # prescribed displacements

# -------------------------
# MATRIX PARTITION
# -------------------------
K_UU = K[0:DOFs, 0:DOFs]
K_UP = K[0:DOFs, DOFs:DOFs+DOCs]
K_PU = K[DOFs:DOFs+DOCs, 0:DOFs]
K_PP = K[DOFs:DOFs+DOCs, DOFs:DOFs+DOCs]

# -------------------------
# SOLVE SYSTEM
# -------------------------
F = Fp - K_UP @ Up
Uu = np.linalg.solve(K_UU, F)

# -------------------------
# REACTION FORCES
# -------------------------
Fu_reactions = K_PU @ Uu + K_PP @ Up

# -------------------------
# UPDATE ENL
# -------------------------
ENL = update_nodes(ENL, Uu, NL, Fu_reactions)

# -------------------------
# OUTPUT (OPTIONAL DEBUG)
# -------------------------
print("Displacements (free DOFs):\n", Uu)
print("Reaction forces:\n", Fu_reactions)











