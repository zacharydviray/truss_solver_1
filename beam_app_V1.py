import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
import tempfile, os, shutil


# -----------------------------
# FEM Beam Element (Euler-Bernoulli)
# -----------------------------
# DOF per node: [w, theta]

class BeamFEM:
    def __init__(self, L, E, I, n_elems=50):
        self.L = L
        self.E = E
        self.I = I
        self.n = n_elems
        self.nnodes = n_elems + 1
        self.dof = 2 * self.nnodes

        self.x = np.linspace(0, L, self.nnodes)
        self.K = np.zeros((self.dof, self.dof))
        self.F = np.zeros(self.dof)

    def beam_k(self, L_e):
        EI = self.E * self.I
        return (EI / L_e**3) * np.array([
            [12, 6*L_e, -12, 6*L_e],
            [6*L_e, 4*L_e**2, -6*L_e, 2*L_e**2],
            [-12, -6*L_e, 12, -6*L_e],
            [6*L_e, 2*L_e**2, -6*L_e, 4*L_e**2]
        ])

    def assemble(self, w_udl=0.0, point_loads=None):
        if point_loads is None:
            point_loads = []

        L_e = self.L / self.n
        k_e = self.beam_k(L_e)

        for e in range(self.n):
            dof_map = [2*e, 2*e+1, 2*e+2, 2*e+3]
            for i in range(4):
                for j in range(4):
                    self.K[dof_map[i], dof_map[j]] += k_e[i, j]

            # consistent UDL load
            if w_udl != 0:
                f_e = np.array([
                    w_udl*L_e/2,
                    w_udl*L_e**2/12,
                    w_udl*L_e/2,
                    -w_udl*L_e**2/12
                ])
                for i in range(4):
                    self.F[dof_map[i]] += f_e[i]

        # point loads
        for P, a in point_loads:
            idx = int(a / self.L * self.n)
            idx = min(max(idx, 0), self.nnodes-1)
            self.F[2*idx] -= P

    def apply_bc(self):
        fixed = [0, 2*(self.nnodes-1)]  # simply supported ends

        for dof in fixed:
            self.K[dof, :] = 0
            self.K[:, dof] = 0
            self.K[dof, dof] = 1
            self.F[dof] = 0

    def solve(self):
        self.U = np.linalg.solve(self.K, self.F)
        self.w = self.U[0::2]
        self.theta = self.U[1::2]
        return self.w

    def reactions(self):
        R = self.K @ self.U - self.F
        return R[0], R[-2]


# -----------------------------
# GUI APP
# -----------------------------

class BeamApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FEM Beam Solver (Improved)")
        self.geometry("850x650")

        self.loads = []
        self.create_widgets()

    def create_widgets(self):
        frame = ttk.LabelFrame(self, text="Inputs")
        frame.pack(fill="x", padx=10, pady=10)

        self.L = self.entry(frame, "Span L (m)", 0)
        self.w = self.entry(frame, "UDL w (kN/m)", 1)
        self.E = self.entry(frame, "E (kN/m²)", 2)
        self.I = self.entry(frame, "I (m⁴)", 3)

        load_frame = ttk.LabelFrame(self, text="Point Loads")
        load_frame.pack(fill="x", padx=10, pady=10)

        self.P = ttk.Entry(load_frame)
        self.a = ttk.Entry(load_frame)

        ttk.Label(load_frame, text="P").grid(row=0, column=0)
        self.P.grid(row=0, column=1)
        ttk.Label(load_frame, text="a").grid(row=0, column=2)
        self.a.grid(row=0, column=3)

        ttk.Button(load_frame, text="Add", command=self.add_load).grid(row=0, column=4)

        self.listbox = tk.Listbox(load_frame, height=5)
        self.listbox.grid(row=1, column=0, columnspan=5, sticky="ew")

        ttk.Button(self, text="Run FEM Analysis", command=self.run).pack(pady=10)

        self.out = tk.Text(self, height=10)
        self.out.pack(fill="both", padx=10, pady=10)

    def entry(self, parent, text, row):
        ttk.Label(parent, text=text).grid(row=row, column=0)
        e = ttk.Entry(parent)
        e.grid(row=row, column=1)
        return e

    def add_load(self):
        try:
            P = float(self.P.get())
            a = float(self.a.get())
            self.loads.append((P, a))
            self.listbox.insert(tk.END, f"P={P} at x={a}")
        except:
            messagebox.showerror("Error", "Invalid load")

    def run(self):
        try:
            L = float(self.L.get())
            w = float(self.w.get())
            E = float(self.E.get())
            I = float(self.I.get())

            print("\n==============================")
            print("RUNNING FEM ANALYSIS")
            print(f"L = {L}, w = {w}, E = {E}, I = {I}")
            print("==============================\n")

            fem = BeamFEM(L, E, I, n_elems=60)
            fem.assemble(w, self.loads)
            fem.apply_bc()
            defl = fem.solve()

            R1, R2 = fem.reactions()

            # LIVE TERMINAL OUTPUT
            print("--- REACTIONS ---")
            print(f"Reaction at node 1 (RA): {R1:.4f} kN")
            print(f"Reaction at node 2 (RB): {R2:.4f} kN")

            print("\n--- SAMPLE INTERNAL RESULTS ---")
            for i in range(0, len(fem.x), max(1, len(fem.x)//10)):
                print(f"x={fem.x[i]:.3f} m | w={defl[i]:.6f} m")

            self.out.delete("1.0", tk.END)
            self.out.insert(tk.END, f"Reaction 1: {R1:.2f}\n")
            self.out.insert(tk.END, f"Reaction 2: {R2:.2f}\n")

            self.plot(fem.x, defl)

        except Exception as e:
            messagebox.showerror("Error", str(e))
        try:
            L = float(self.L.get())
            w = float(self.w.get())
            E = float(self.E.get())
            I = float(self.I.get())

            fem = BeamFEM(L, E, I, n_elems=60)
            fem.assemble(w, self.loads)
            fem.apply_bc()
            defl = fem.solve()

            R1, R2 = fem.reactions()

            self.out.delete("1.0", tk.END)
            self.out.insert(tk.END, f"Reaction 1: {R1:.2f}\n")
            self.out.insert(tk.END, f"Reaction 2: {R2:.2f}\n")

            self.plot(fem.x, defl)

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def plot(self, x, w):
        plt.figure()
        plt.plot(x, w)
        plt.title("Deflection (FEM)")
        plt.xlabel("x")
        plt.ylabel("w")
        plt.grid()
        plt.show()


if __name__ == "__main__":
    app = BeamApp()
    app.mainloop()
