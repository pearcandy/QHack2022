import sys
import pennylane as qml
from pennylane import numpy as np
from pennylane import hf


def ground_state_VQE(H):
    """Perform VQE to find the ground state of the H2 Hamiltonian.

    Args:
        - H (qml.Hamiltonian): The Hydrogen (H2) Hamiltonian

    Returns:
        - (float): The ground state energy
        - (np.ndarray): The ground state calculated through your optimization routine
    """

    # QHACK #
    qubits = 4
    dev = qml.device("default.qubit", wires=qubits)

    def circuit(param, wires):
        qml.BasisState(np.array([1, 1, 0, 0]), wires=wires)
        qml.SingleExcitation(param[0], wires=[0, 2])
        qml.SingleExcitation(param[1], wires=[1, 3])
        qml.DoubleExcitation(param[2], wires=[0, 1, 2, 3])

    @qml.qnode(dev)
    def cost_fn(param):
        circuit(param, wires=range(qubits))
        return qml.expval(H)

    opt = qml.GradientDescentOptimizer(stepsize=0.4)
    theta = np.array([np.random.random() for i in range(3)],
                     requires_grad=True)

    # store the values of the cost function
    energy = [cost_fn(theta)]

    # store the values of the circuit parameter
    angle = [theta]

    max_iterations = 100
    conv_tol = 1e-07

    for n in range(max_iterations):
        theta, prev_energy = opt.step_and_cost(cost_fn, theta)

        energy.append(cost_fn(theta))
        angle.append(theta)

        conv = np.abs(energy[-1] - prev_energy)

        if n % 2 == 0:
            print(f"Step = {n},  Energy = {energy[-1]:.8f} Ha")

        if conv <= conv_tol:
            break

    #|1100>
    theta1 = angle[-1][0]
    theta2 = angle[-1][1]
    theta3 = angle[-1][2]

    g3 = np.cos(theta1 / 2.0) * np.cos(theta2 / 2.0) * np.cos(
        theta3 / 2.0) + np.sin(theta1 / 2.0) * np.sin(theta2 / 2.0) * np.sin(
            theta3 / 2.0)

    g6 = -np.sin(theta1 / 2.0) * np.cos(theta2 / 2.0)
    g9 = -np.cos(theta1 / 2.0) * np.sin(theta2 / 2.0)
    g12 = np.sin(theta1 / 2.0) * np.sin(theta2 / 2.0) * np.cos(
        theta3 / 2.0) - np.cos(theta1 / 2.0) * np.cos(theta2 / 2.0) * np.sin(
            theta3 / 2.0)

    ground_state = np.array(
        [0, 0, 0, g3, 0, 0, g6, 0, 0, g9, 0, 0, g12, 0, 0, 0])
    return energy[-1], ground_state

    # QHACK #


def create_H1(ground_state, beta, H):
    """Create the H1 matrix, then use `qml.Hermitian(matrix)` to return an observable-form of H1.

    Args:
        - ground_state (np.ndarray): from the ground state VQE calculation
        - beta (float): the prefactor for the ground state projector term
        - H (qml.Hamiltonian): the result of hf.generate_hamiltonian(mol)()

    Returns:
        - (qml.Observable): The result of qml.Hermitian(H1_matrix)
    """

    # QHACK #
    H0 = np.outer(ground_state, np.conj(ground_state))
    #Hmat = qml.utils.sparse_hamiltonian(H).real.toarray()
    Hmat = H.matrix
    H1_mat = Hmat + beta * H0
    return qml.Hermitian(H1_mat, wires=[0, 1, 2, 3])
    # QHACK #


def excited_state_VQE(H1):
    """Perform VQE using the "excited state" Hamiltonian.

    Args:
        - H1 (qml.Observable): result of create_H1

    Returns:
        - (float): The excited state energy
    """

    # QHACK #
    qubits = 4
    dev = qml.device("default.qubit", wires=qubits)

    def circuit2(param, wires):
        qml.BasisState(np.array([1, 1, 0, 0]), wires=wires)
        qml.SingleExcitation(param[0], wires=[0, 2])
        qml.SingleExcitation(param[1], wires=[1, 3])
        qml.DoubleExcitation(param[2], wires=[0, 1, 2, 3])

    @qml.qnode(dev)
    def cost_fn2(param):
        circuit2(param, wires=range(qubits))
        return qml.expval(H1)

    opt = qml.GradientDescentOptimizer(stepsize=0.2)
    theta = np.array([np.random.random() for i in range(3)],
                     requires_grad=True)

    # store the values of the cost function
    energy = [cost_fn2(theta)]

    # store the values of the circuit parameter
    angle = [theta]

    max_iterations = 100
    conv_tol = 1e-07

    for n in range(max_iterations):
        theta, prev_energy = opt.step_and_cost(cost_fn2, theta)

        energy.append(cost_fn2(theta))
        angle.append(theta)

        conv = np.abs(energy[-1] - prev_energy)

        if n % 2 == 0:
            print(f"Step = {n},  Energy = {energy[-1]:.8f} Ha")

        if conv <= conv_tol:
            break

    print(energy[-1], angle[-1])
    return energy[-1]

    # QHACK #


if __name__ == "__main__":
    coord = float(sys.stdin.read())
    symbols = ["H", "H"]
    geometry = np.array([[0.0, 0.0, -coord], [0.0, 0.0, coord]],
                        requires_grad=False)
    mol = hf.Molecule(symbols, geometry)

    H = hf.generate_hamiltonian(mol)()
    E0, ground_state = ground_state_VQE(H)

    #beta = 0.0
    beta = 15.0
    H1 = create_H1(ground_state, beta, H)
    E1 = excited_state_VQE(H1)

    answer = [np.real(E0), E1]
    print(*answer, sep=",")
