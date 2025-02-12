import numpy as np
from scipy.linalg import expm

Z0 = 377 #impedance of free space (Ohms)
class Berreman(object):
    def __init__(self):
        pass

    def rxmat(self, x):
        """
            Rotation matrix around the x-axis.

            Parameters:
                x: Rotation angle in radians.

            Returns:
                y: A 3x3 rotation matrix for rotation around the x-axis.
            """
        c = np.cos(x)
        s = np.sin(x)
        y = np.array([
            [1, 0, 0],
            [0, c, -s],
            [0, s, c]
        ])
        return y

    def rzmat(self, x):
        """
        Rotation matrix around the z-axis.
        Parameters:
            x: Rotation angle in radians.
        Returns:
            y: A 3x3 rotation matrix for rotation around the z-axis.
        """
        c = np.cos(x)
        s = np.sin(x)
        y = np.array([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])
        return y

    def rbmat(self,xi):
        """
        Rotation matrix used in Berreman calculus.
        :param xi:
        :return: the rotation matix for angle xi radian
        """
        c = np.cos(xi)
        s = np.sin(xi)
        Rb = np.array([[c, 0, -s, 0],
              [0, c, 0, s],
              [s, 0, c, 0],
              [0, -s, 0, c]])
        return Rb

    def epsilon(self, layer):
        """
        Relative permittivity matrix used in Berreman calculus.
        Parameters:
        layer: Array-like of size (1, 6), containing:
               [n1, n2, n3, eta, psi, xi], where:
               - n1, n2, n3: Principal refractive indices.
               - eta, psi, xi: Euler angles in radians.
        Returns:
        y: The 3x3 relative permittivity matrix.
        """
        if len(layer) == 6:
            # Diagonal matrix with squared refractive indices
            e = np.diag([layer[0] ** 2, layer[1] ** 2, layer[2] ** 2])

            # Rotate around the x-axis
            e = self.rxmat(layer[3]) @ e @ self.rxmat(-layer[3])

            # Rotate around the z-axis
            e = self.rzmat(layer[4]) @ e @ self.rzmat(-layer[4])

            # Rotate around the x-axis again
            y = self.rxmat(layer[5]) @ e @ self.rxmat(-layer[5])

            return y
        else:
            raise ValueError("Input layer must have 6 elements.")

    def fmat(self, material, beta):
        """
        Field matrix, alphas and fields used in Berreman calculus.

        material: list describing the material properties.
                                  - 1 element: Isotropic material (e.g., refractive index).
                                  - 2 elements: PS medium (orthogonal polarizations).
                                  - 6 elements: Birefringent anisotropic material (detailed parameters).
        beta (float): Propagation parameter. Use 0 for PS medium.
                    Otherwise, Propagation constant
                    k0 = 2 * pi / lambda; % Wavenumber in vacuum
                    beta = k0 * n * sind(theta_i);

        :return:
        F (numpy.ndarray): Field matrix.
        alpha (numpy.ndarray): Propagation constants.
        E (numpy.ndarray): Electric field matrix.
        H (numpy.ndarray): Magnetic field matrix.
        """
        if not isinstance(material, list): # Isotropic material
            n = material
            alpha_val = np.sqrt(n**2 - beta**2) # alpha = ncos(theta)
            gp = n ** 2 / (alpha_val * Z0)
            gs = -alpha_val / Z0
            F = np.array([[1, 1, 0, 0],
                          [gp, -gp, 0, 0],
                          [0, 0, 1, 1],
                          [0, 0, gs, -gs]])

            alpha = np.array([alpha_val, -alpha_val, alpha_val, -alpha_val])

            E = np.array([-Z0 * beta * F[1] / (n ** 2),
                         F[0],
                         F[2]])

            H = np.array([beta * F[2] / Z0,
                         F[3],
                         F[1]])

        elif len(material) == 2 and beta == 0:  # PS material
            print("PS material")
            ny = material[0]
            nz = material[1]

            F = np.array([[1, 1, 0, 0],
                          [ny / Z0, -ny / Z0, 0, 0],
                          [0, 0, 1, 1],
                          [0, 0, -nz / Z0, nz / Z0]])
            alpha = np.array([ny, -ny, nz, -nz])
            E = np.array([[0, 0, 0, 0],
                          F[0],
                          F[2]])
            H = np.array([[0, 0, 0, 0],
                          F[3],
                          F[1]])

        elif len(material) == 6:  # Birefringent material
            n1, n2, n3, eta, psi, xi = material
            if abs(n1 - n2) < 1e-6 and abs(n1 - n3) < 1e-6:  # Test for isotropic layer
                alpha_val = np.sqrt(n1 ** 2 - beta ** 2)
                gs = -alpha_val / Z0
                gp = n1 ** 2 / (alpha_val * Z0)
                F = np.array([[1, 1, 0, 0],
                              [gp, -gp, 0, 0],
                              [0, 0, 1, 1],
                              [0, 0, gs, -gs]])
                alpha = np.array([alpha_val, -alpha_val, alpha_val, -alpha_val])
                E = np.array([-Z0 * beta * F[1] / (n1 ** 2),
                              F[0],
                              F[2]])
                H = np.array([beta * F[2] / Z0,
                              F[3],
                              F[1]])
            else:
                e = self.epsilon(material)  # Dielectric matrix function
                exx, exy, exz = e[0]
                eyy, eyz = e[1, 1], e[1, 2]
                ezz = e[2, 2]

                Mbeta = np.array([
                    [-beta * exy / exx, Z0 * (1 - beta ** 2 / exx), -beta * exz / exx, 0],
                    [(eyy - exy ** 2 / exx) / Z0, -beta * exy / exx, (eyz - exy * exz / exx) / Z0, 0],
                    [0, 0, 0, -Z0],
                    [(-eyz + exz * exy / exx) / Z0, beta * exz / exx, (beta ** 2 - ezz + exz ** 2 / exx) / Z0, 0]
                ])
                # Compute eigenvalues and eigenvectors:
                eigenvalues, eigenvectors = np.linalg.eig(Mbeta)
                #Ensures eigenvector is 2D
                alpha = eigenvalues #diagonal matrix

                # Normalize columns of F
                normc = np.sqrt(np.abs(eigenvectors[0] ** 2) + np.abs(eigenvectors[2] ** 2)) * np.exp(1j * np.angle(eigenvectors[0] + eigenvectors[2]))
                F = np.real(np.array(eigenvectors)/normc)

                # Sort columns of F based on alpha #TODO this is working but still don't fully understand why
                sorted_indices = np.argsort(np.abs(alpha))
                alpha = alpha[-sorted_indices + 1]
                F = F[:, -sorted_indices + 1]

                alpha = np.around(alpha, 4)
                F = np.around(F, 4)

                E = np.array([-(exy * F[0] + exz * F[2] + Z0 * beta * F[1]) / exx,
                               F[0],
                               F[2]])

                E = np.around(E, 4) #todo: slight different values in the first row of the example, but should be OK

                H = np.array([beta * F[2] / Z0, F[3], F[1]])

                H = np.round(H, 4)

        else:
            raise ValueError("Invalid material matrix or beta parameter.")

        return [F, alpha, E, H]

    def pmat(self, alpha, dw):
        """
        Phase matrix used in Berreman calculus.
        Parameters:
            alpha: Array-like of propagation constants (length 4).
            dw: Layer thickness.
        Returns:
        Ad: Phase matrix (4x4).
        """
        c = -1j * 2 * np.pi * dw
        Ad = np.diag([
            np.exp(c * alpha[0]),
            np.exp(c * alpha[1]),
            np.exp(c * alpha[2]),
            np.exp(c * alpha[3])
        ])
        return Ad

    def cmat(self, arg1, beta):
        """
        Characteristic matrix in Berreman calculus.
        Supports various formats:
        - cmat(layer, beta)
        - cmat(stack, beta)
        - cmat(system, beta)
        - cmat([layer, gain], beta) for a PS medium with single-pass amplitude gain.

    The input parameter beta must be entered as 0 for a PS medium
        :param arg1: Input matrix describing the system, stack, or layer.
        :param beta: Propagation parameter (0 for a PS medium).
        :return: M: Characteristic matrix (4x4).
        """

        arg1 = np.array(arg1)

        if arg1.ndim == 1:
            arg1 = arg1.reshape(1, -1)

        rows, cols = arg1.shape
        M = np.eye(4)
        if cols == 7:  # Birefringent case
            m1, m2 = 1, rows
            if np.isinf(arg1[0, 6]) and np.isinf(arg1[-1, 6]):
                m1, m2 = 2, rows - 1

            for j in range(m1 - 1, m2):  # Adjusted for 0-based indexing in Python
                F, alpha, E, H = self.fmat(list(arg1[j][:6]), beta)
                M = M @ F @ self.pmat(alpha, arg1[j, 6]) @ np.linalg.inv(F)

        elif cols == 2:  # Isotropic case
            m1, m2 = 1, rows
            if np.isinf(arg1[0, 1]) and np.isinf(arg1[-1, 1]):
                m1, m2 = 2, rows - 1
            for j in range(m1 - 1, m2):
                n = arg1[j, 0]
                alph = np.sqrt(n**2 - beta**2)
                gp = n * n / alph / Z0
                gs = -alph / Z0
                dw = arg1[j, 1]
                phi = 2 * np.pi * dw * alph
                c = np.cos(phi)
                s = np.sin(phi)
                M = M @ np.array([
                    [c, -1j / gp * s, 0, 0],
                    [-1j * gp * s, c, 0, 0],
                    [0, 0, c, -1j / gs * s],
                    [0, 0, -1j * gs * s, c]
                ])

        elif cols == 3 and beta == 0:  # PS case
            m1, m2 = 1, rows
            if np.isinf(arg1[0, 2]) and np.isinf(arg1[-1, 2]):
                m1, m2 = 2, rows - 1
            for j in range(m1 - 1, m2):
                ny = arg1[j, 0]
                nz = arg1[j, 1]
                dw = arg1[j, 2]
                phiy = 2 * np.pi * ny * dw
                phiz = 2 * np.pi * nz * dw
                cy, cz = np.cos(phiy), np.cos(phiz)
                sy, sz = np.sin(phiy), np.sin(phiz)
                M = M @ np.array([
                    [cy, -1j / ny * Z0 * sy, 0, 0],
                    [-1j * ny / Z0 * sy, cy, 0, 0],
                    [0, 0, cz, 1j / nz * Z0 * sz],
                    [0, 0, 1j * nz / Z0 * sz, cz]
                ])
        elif rows == 1 and cols == 4 and beta == 0:  # PS layer with gain
            ny, nz, dw, g = arg1[0, :]
            phiy = 2 * np.pi * ny * dw
            phiz = 2 * np.pi * nz * dw
            F = np.array([
                [1, 1, 0, 0],
                [ny / Z0, -ny / Z0, 0, 0],
                [0, 0, 1, 1],
                [0, 0, -nz / Z0, nz / Z0]
            ])
            Ad = np.diag([
                np.exp(-1j * phiy) / g,
                np.exp(1j * phiy) * g,
                np.exp(-1j * phiz) / g,
                np.exp(1j * phiz) * g
            ])
            M = F @ Ad @ np.linalg.inv(F)

        else:
            raise ValueError("Unsupported input format for cmat.")

        return M

    def cmmat(self, Fc, M, Fs):
        """
        Mueller and Stokes matrices for a birefringent coating, including the effects
        of interference.
        :param Fc: Field Matrix of the cover
        :param M:
        :param Fs: Field Matrix of the substrate
        :return:
            Mr: 4x4 Mueller matrix of the reflected beam
            Mt: 4x4 Mueller matrix of the transmitted beam
            sR: 4x6 Stokes matrix of the reflected beam
            sT: 4x6 Stokes matrix of the transmitted beam
        """
        # Equation 5.31 synthetic matrix of stokes vectors
        S = np.array([
            [1, 1, 1, 1, 1, 1],
            [1, -1, 0, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 0, 1, -1]
        ])

        gCp, gSp = Fc[0][1, 0], Fs[0][1, 0]
        gCs, gSs = Fc[0][3, 2], Fs[0][3, 2]

        ctC = np.sqrt(-gCs / gCp)
        ctS = np.sqrt(-gSs / gSp)

        R, r = self.reflect(Fc, M, Fs, 0, 0)

        r11, r12 = r[0, :2]
        r21, r22 = r[1, :2]
        t11, t12 = r[2, :2]
        t21, t22 = r[3, :2]

        ## Computing the different polarization states
        # POLARIZATION STATES DEFINITION:
        POLAR = [[1, 0],  # P0
                 [0, 1],  # P90
                 [1 / np.sqrt(2), 1 / np.sqrt(2)],  # P45
                 [1 / np.sqrt(2), - 1 / np.sqrt(2)],  # P135
                 [1 / np.sqrt(2), -1j / np.sqrt(2)],  # R
                 [1 / np.sqrt(2), +1j / np.sqrt(2)]  # L
                 ]
        sr_, st_ = [], []  # stores the stokes components for the reflected and transmited beams

        # TODO: works for DKIST example but not with the book example
        for state in POLAR:
            EY, EZ = state
            Ey = EY * ctC
            Ez = EZ
            EyC = np.outer(r11, Ey) + np.outer(r12, Ez)
            EzC = np.outer(r21, Ey) + np.outer(r22, Ez)
            EY = -EyC / ctC
            EZ = EzC
            sr_.append(self.stokes(EY, EZ))

            EyS = np.outer(t11, Ey) + np.outer(t12, Ez)
            EzS = np.outer(t21, Ey) + np.outer(t22, Ez)
            EY = EyS / ctS
            EZ = EzS
            st = self.stokes(EY, EZ)
            st = np.outer(st, (gSs / gCs))
            st_.append(st)

        # Builds the transposed array of 4x6 stokes Matrix for the reflected and transmitted beam
        sR = np.array(sr_)
        sT = np.array(st_)
        sR = sR[:, :, 0, 0].T
        sT = sT[:, :, 0].T
        S_inv = np.linalg.pinv(S)

        # Equation 5.34, builds the 4x4 Mueller Matrix for the reflected and transmitted beam
        Mr = np.dot(sR, S_inv)
        Mt = np.dot(sT, S_inv)

        return Mr, Mt, sR, sT

    def hmat(self, arg1, beta=None):
        """
        Characteristic matrix of a chiral film for use in Berreman calculus.
        :return:
        """
        cols = len(arg1)

        if beta is None and cols == 6:  # Case 1: Double helix, PS sublayers
            n2, n3, xi1, xi2, nj, G = arg1
            nav = (abs(n2) + abs(n3)) / 2
            dxi = xi2 - xi1

            if nj == 0:  # Continuous matrix
                M = expm(np.array([
                    [0, -1j * G * abs(dxi) / nav, -dxi, 0],
                    [-1j * G * abs(dxi) * n2 ** 2 / nav, 0, 0, dxi],
                    [dxi, 0, 0, 1j * G * abs(dxi) / nav],
                    [0, -dxi, 1j * G * abs(dxi) * n3 ** 2 / nav, 0]
                ]))
            else:  # Deposition matrix
                dxij = dxi / nj
                dw = G / (2 * nav * nj) * abs(dxi) / np.pi
                M = self.cmat([n2, n3, dw / 2], 0)
                M = np.linalg.matrix_power(np.matmul(np.matmul(M, self.rbmat(dxij)), M), nj)

            M = np.matmul(np.matmul(self.rbmat(xi1), M), self.rbmat(-xi2))

        elif beta is None and cols == 7:  # Case 2: Double helix with gain
            n2, n3, xi1, xi2, nj, G, g = arg1
            nav = (abs(n2) + abs(n3)) / 2
            dxi = (xi2 - xi1) / nj
            dw = G / (2 * nav * nj) * abs(xi2 - xi1) / np.pi

            M = self.cmat([n2, n3, dw / 2, g], 0)
            M = np.linalg.matrix_power(np.matmul(np.matmul(M, self.rbmat(dxi)), M), nj)
            M = np.matmul(np.matmul(self.rbmat(xi1), M), self.rbmat(-xi2))

        elif beta is not None and cols == 9:  # Case 3: Single helix, tilted sublayers
            n1, n2, n3, eta, psi, xi1, xi2, nj, G = arg1
            np_val = 1 / np.sqrt((np.sin(psi) ** 2) / (n1 ** 2) + (np.cos(psi) ** 2) / (n2 ** 2))
            nav = (abs(np_val) + abs(n3)) / 2
            dxi = (xi2 - xi1) / nj
            dw = G / (2 * nav * nj) * abs(xi2 - xi1) / np.pi

            if beta == 0:
                M = self.cmat([n1, n2, n3, eta, psi, 0, dw / 2], 0)
                M = np.linalg.matrix_power(np.matmul(np.matmul(M, self.rbmat(dxi)), M), nj)
                M = np.matmul(np.matmul(self.rbmat(xi1), M), self.rbmat(-xi2))

            else:
                M = self.cmat([n1, n2, n3, eta, psi, xi1, dw / 2], beta)
                for j in range(2, nj + 1):
                    xi = xi1 + (xi2 - xi1) * (j - 1) / nj
                    M = np.matmul(M, self.cmat([n1, n2, n3, eta, psi, xi, dw], beta))
                M = np.matmul(M, self.cmat([n1, n2, n3, eta, psi, xi2, dw / 2], beta))
        else:
            raise ValueError("Invalid input arguments.")
        return M

    def poynting(self, E, H = None):
        """
        Poynting flux vectors used in Berreman calculus.
        Allowed formats are poynting(E, H) and poynting(F).

        Parameters:
            E: Electric field matrix (numpy.ndarray).
            H: Magnetic field matrix (numpy.ndarray), optional.

        Returns:
            p: 3x4 pointing Matrix in which the rows correspond to px, py and pz
        """
        if H is not None:  # When both E and H are provided corresponds to Eq 3.55
            p = np.array([
                np.real(E[1] * np.conj(H[2]) - E[2] * np.conj(H[1])),
                np.real(E[2] * np.conj(H[0]) - E[0] * np.conj(H[2])),
                np.real(E[0] * np.conj(H[1]) - E[1] * np.conj(H[0]))
            ]) / 2
        else:  # When only F (which is equivalent to E in this case) is provided
            F = E
            p = np.real((F[0] * np.conj(F[1]) - F[2] * np.conj(F[3]))) / 2

            #Funciona!
        return np.round(p, 4)

    def smat(self, arg1, arg2, arg3):
        """
            System matrix used in Berreman calculus.
            Formats are smat(Fc, M, Fs) and smat(system, beta).
            The input parameter beta must be entered as 0 for a PS medium.

            Parameters:
                arg1: Either the field matrix `Fc` (numpy.ndarray), or `system` (list or numpy.ndarray).
                arg2: If `arg1` is a `system`, this is the `beta` (float).
                arg3: If provided, this is the `Fs` (numpy.ndarray).
            Returns:
                A: System matrix result (numpy.ndarray).
            """
        if arg3 is not None:  # smat(Fc, M, Fs)
            Fc, M, Fs = arg1, arg2, arg3

        elif len(arg1.shape) == 2 and (
                arg1.shape[1] == 2 or arg1.shape[1] == 7):  # smat(system, beta) birefringent or isotropic
            system, beta = arg1, arg2
            rows, cols = system.shape
            Fc = self.fmat(system[0, :cols - 1], beta)
            M = self.cmat(system, beta)
            Fs = self.fmat(system[rows - 1, :cols - 1], beta)

        elif len(arg1.shape) == 2 and arg1.shape[1] == 3 and arg2 == 0:  # smat(system, beta) PS medium
            system, beta = arg1, arg2
            rows, cols = system.shape
            Fc = self.fmat(system[0, :cols - 1], beta)
            M = self.cmat(system, beta)
            Fs = self.fmat(system[rows - 1, :cols - 1], beta)

        else:
            raise ValueError("Invalid arguments")
        # Compute the system matrix result
        Fc_inv = np.linalg.inv(Fc[0])
        A = (Fc_inv @ M) @ Fs[0]
        return np.round(A, 5)

    def stokes(self, Ey, Ez):
        """
                    Calculate the Stokes parameters from the electric field components Ey and Ez.

                    Parameters:
                        Ey (complex): Ey component of the electric field.
                        Ez (complex): Ez component of the electric field.

                    Returns:
                        tuple: A tuple containing Stokes parameters (S0, S1, S2, S3).
                    """
        # Compute the phase difference
        Del = np.angle(Ey) - np.angle(Ez)

        # Stokes parameters
        S0 = abs(Ey) ** 2 + abs(Ez) ** 2
        S1 = abs(Ey) ** 2 - abs(Ez) ** 2
        S2 = 2 * abs(Ey) * abs(Ez) * np.cos(Del)
        S3 = 2 * abs(Ey) * abs(Ez) * np.sin(Del)

        return S0, S1, S2, S3

    def reflect(self, *args):
        """
    Reflection and transmission coefficients in Berreman calculus.
    Supports multiple input formats as described below:

    Usage:
        R, r = reflect(Fc, M, Fs)
        R, r = reflect(Fc, M, Fs, e, phi)
        R, r = reflect(system, beta)
        R, r = reflect(system, beta, e, phi)
    Parameters:
        args: Variable-length arguments based on usage.
              - `Fc`, `M`, `Fs`: Field matrices for input/output and system.
              - `system`: A matrix describing material layers.
              - `beta`: Propagation parameter.
              - `e`, `phi`: Ellipticity and azimuthal angle (optional).
        :return:
                R: Reflection matrix.
                r: Transmission matrix.
        """
        I = np.eye(4) # 4x4 identity matrix

        ms = len(args[0])

        if len(args) == 3 and ms == 4:  # reflect(Fc, M, Fs)
            Fc, M, Fs = args
            e, phi = 0, 0
            print("reflect(Fc, M, Fs)")

        elif len(args) == 5 and ms == 4:  # reflect(Fc, M, Fs, e, phi)
            Fc, M, Fs, e, phi = args
            nC = np.sqrt(-Fc[0][1][0] * Fc[0][3][2] * Z0 * Z0)
            beta = nC * np.sqrt(1 - nC ** 2 / (Fc[0][1][0] ** 2 * Z0 ** 2))
            ctC = np.sqrt(1 - beta ** 2 / nC ** 2)
            nS = np.sqrt(-Fs[0][1][0] * Fs[0][3][2] * Z0 ** 2)
            ctS = np.sqrt(1 - beta ** 2 / nS ** 2)
            print("reflect(Fc, M, Fs, e, phi)")

        elif len(args) == 2 and (ms == 2 or ms == 7):  # reflect(system, beta) birefringent or isotropic
            system, beta = args
            rows, cols = np.asarray(system).shape
            Fc = self.fmat(system[0, :cols - 1], beta)
            M = self.cmat(system, beta)
            Fs = self.fmat(system[rows - 1, :cols - 1], beta)
            print("reflect(system, beta) birefringent or isotropic")

        elif len(args) == 2 and ms == 3 and args[1] == 0:  # reflect(system, beta) PS
            system, beta = args
            rows, cols = system.shape
            Fc = self.fmat(system[0, :cols - 1], beta)
            M = self.cmat(system, beta)
            Fs = self.fmat(system[rows - 1, :cols - 1], beta)
            print("# reflect(system, beta) PS")

        elif len(args) == 4 and (ms == 2 or ms == 7):  # reflect(system, beta, e, phi) birefringent or isotropic
            system, beta, e, phi = args
            rows, cols = system.shape
            Fc = self.fmat(system[0, :cols - 1], beta)
            M = self.cmat(system, beta)
            Fs = self.fmat(system[rows - 1, :cols - 1], beta)
            nC = system[0, 0]
            nS = system[rows - 1, 0]
            ctC = np.sqrt(1 - beta ** 2 / nC ** 2)
            ctS = np.sqrt(1 - beta ** 2 / nS ** 2)
            print("reflect(system, beta, e, phi) birefringent or isotropic")

        elif len(args) == 4 and ms == 3 and args[1] == 0:  # reflect(system, beta, e, phi) PS
            system, beta, e, phi = args
            rows, cols = system.shape
            Fc = self.fmat(system[0, :cols - 1], beta)
            M = self.cmat(system, beta)
            Fs = self.fmat(system[rows - 1, :cols - 1], beta)
            nC = system[0, 0]
            nS = system[rows - 1, 0]
            ctC = np.sqrt(1 - beta ** 2 / nC ** 2)
            ctS = np.sqrt(1 - beta ** 2 / nS ** 2)
            print("# reflect(system, beta, e, phi) PS")

        else:
            raise ValueError("Invalid combination of inputs for reflect function.")

        # p’s are the powers carried by the basis vectors along the x-axis, which is the first row in case a 4x4 matrix
        # is returned
        pc = self.poynting(Fc[0])
        ps = self.poynting(Fs[0])
        cp = 0.00001

        for i in [0, 2]:
            if abs(pc[i]) < cp:
                pc[i] = np.nan
        for i in [1, 3]:
           if abs(pc[i]) < cp:
               ps[i] = np.nan

        ## PUT CHANGE BASIS VECTOR HERE, DEPRECIATED TODO: COMPROBAR CODIGO ###
        if e != 0 or phi != 0:
            c = np.cos(phi)
            s = np.sin(phi)
            print("YASSSSS")

            # Intermediate variables
            gCp = nC / (ctC * Z0)
            gCs = -nC * ctC / Z0

            # Field matrix for `Fc`
            Fc = np.array([
                [(c + 1j * e * s) * ctC, (c - 1j * e * s) * ctC, -(s + 1j * e * c) * ctC, -(s - 1j * e * c) * ctC],
                [gCp * (c + 1j * e * s) * ctC, -gCp * (c - 1j * e * s) * ctC, -gCp * (s + 1j * e * c) * ctC,
                 gCp * (s - 1j * e * c) * ctC],
                [(s - 1j * e * c), (s + 1j * e * c), (c - 1j * e * s), (c + 1j * e * s)],
                [gCs * (s - 1j * e * c), -gCs * (s + 1j * e * c), gCs * (c - 1j * e * s), -gCs * (c + 1j * e * s)]
            ])

            # Intermediate variables
            gSp = nS / (ctS * Z0)
            gSs = -nS * ctS / Z0

            # Field matrix for `Fs`
            Fs = np.array([
                [(c + 1j * e * s) * ctS, (c - 1j * e * s) * ctS, -(s + 1j * e * c) * ctS, -(s - 1j * e * c) * ctS],
                [gSp * (c + 1j * e * s) * ctS, -gSp * (c - 1j * e * s) * ctS, -gSp * (s + 1j * e * c) * ctS,
                 gSp * (s - 1j * e * c) * ctS],
                [(s - 1j * e * c), (s + 1j * e * c), (c - 1j * e * s), (c + 1j * e * s)],
                [gSs * (s - 1j * e * c), -gSs * (s + 1j * e * c), gSs * (c - 1j * e * s), -gSs * (c + 1j * e * s)]
            ])
            ###### DEPRECIATED FOR THE MOMENT BECAUSE WE DO NOT USE IT

        A = self.smat(Fc, M, Fs)

        # Constructing the matrix for inversion
        mat = np.column_stack([I[:, 1], I[:, 3], -A[:, 0], -A[:, 2]])
        cond_number = np.linalg.cond(mat)
        if cond_number > 1e12:
            print("Matrix is ill-conditioned; using pseudoinverse.")
            mat_inv = np.linalg.pinv(mat)
        else:
            mat_inv = np.linalg.inv(mat)
        # Constructing the right-hand side
        rhs = np.column_stack([-I[:, 0], -I[:, 2], A[:, 1], A[:, 3]])

        # Solving for r as per Eq 5.6
        r = mat_inv @ rhs

        # Compute weights and R

        # weights = np.array([
        #     [ps[0]/pc[0], ps[0]/pc[2], ps[0]/ps[2], ps[0]/ps[3]],
        #     [ps[1]/pc[0], ps[1]/pc[1], ps[1]/ps[2], ps[1]/ps[3]],
        #     [pc[2]/pc[0], pc[2]/pc[1], pc[2]/ps[2], pc[2]/ps[3]],
        #     [pc[3]/pc[0], pc[3]/pc[1], pc[3]/ps[2], pc[3]/ps[3]]
        # ])
        weights = np.outer(np.array([pc[1], pc[3], ps[0], ps[2]]), np.array([1/pc[0], 1/pc[2], 1/ps[1], 1/ps[3]]))

        R = np.abs(weights) * np.abs(r) ** 2 # Corresponds to Eq. 5.9

        return R, r