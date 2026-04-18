import numpy as np
import matplotlib.pyplot as plt
from truss_solver_functions import *

# -------------------------
# INPUT
# -------------------------

#NL = np.array([[0, 0],
 #              [1, 0],
  #             [0.5, 1]])

#EL = np.array([[1, 2],
 #              [2, 3],
  #             [3, 1]])

#DorN = np.array([[-1, -1],
 #                [1, -1],
  #               [1, 1]])

#F_input = np.array([[0, 0],
 #                   [0, 0],
  #                  [0, -20]])

#U_input = np.array([[0, 0],
 #                   [0, 0],
  #                  [0, 0]])

def read_nodes():
    nodes = []
    while True:
        line = input("Enter node x y (blank to finish): ").strip()
        if not line:
            break
        x, y = map(float, line.split())
        nodes.append([x, y])
    return np.array(nodes, dtype=float)

def read_elements():
    elements = []
    while True:
        line = input("Enter element n1 n2 (blank to finish): ").strip()
        if not line:
            break
        n1, n2 = map(int, line.split())
        elements.append([n1, n2])
    return np.array(elements, dtype=int)

def read_boundary_conditions(NoN):
    print("Enter boundary conditions for each node (1 for free, -1 for constrained):")
    DorN = []
    for i in range(NoN):
        while True:
            try:
                line = input(f"Node {i+1} (x y): ").strip()
                dx, dy = map(int, line.split())
                DorN.append([dx, dy])
                break
            except ValueError:
                print("Invalid input. Enter two integers separated by space.")
    return np.array(DorN, dtype=int)

def read_forces(NoN):
    print("Enter applied forces for each node (Fx Fy):")
    F_input = []
    for i in range(NoN):
        while True:
            try:
                line = input(f"Node {i+1} (Fx Fy): ").strip()
                fx, fy = map(float, line.split())
                F_input.append([fx, fy])
                break
            except ValueError:
                print("Invalid input. Enter two floats separated by space.")
    return np.array(F_input, dtype=float)

def read_displacements(NoN):
    print("Enter prescribed displacements for each node (Ux Uy, 0 if free):")
    U_input = []
    for i in range(NoN):
        while True:
            try:
                line = input(f"Node {i+1} (Ux Uy): ").strip()
                ux, uy = map(float, line.split())
                U_input.append([ux, uy])
                break
            except ValueError:
                print("Invalid input. Enter two floats separated by space.")
    return np.array(U_input, dtype=float)

def plot_truss(NL, EL):
    plt.figure(figsize=(8, 6))
    
    # Plot nodes
    plt.scatter(NL[:, 0], NL[:, 1], c='blue', s=100, zorder=5)
    
    # Label nodes
    for i, (x, y) in enumerate(NL):
        plt.text(x, y + 0.05, f'{i+1}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Plot elements
    for elem in EL:
        n1, n2 = elem - 1  # Adjust for 0-based indexing
        x1, y1 = NL[n1]
        x2, y2 = NL[n2]
        plt.plot([x1, x2], [y1, y2], 'k-', linewidth=2)
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Truss Structure')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()


    

NL = read_nodes()
EL = read_elements()

# Plot the truss structure
plot_truss(NL, EL)

DorN = read_boundary_conditions(len(NL))
F_input = read_forces(len(NL))
U_input = read_displacements(len(NL))

E = 1e6
A = 0.01

PD = NL.shape[1]
NoN = NL.shape[0]

def truss_check(NL, EL)
    


# plots what the truss would look like

# checks boundary conditions

# asks if this is the final truss to analyze.

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











