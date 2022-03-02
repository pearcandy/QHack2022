#! /usr/bin/python3

import sys
import numpy as np


def givens_rotations(a, b, c, d):
    """Calculates the angles needed for a Givens rotation to out put the state with amplitudes a,b,c and d

    Args:
        - a,b,c,d (float): real numbers which represent the amplitude of the relevant basis states (see problem statement). Assume they are normalized.

    Returns:
        - (list(float)): a list of real numbers ranging in the intervals provided in the challenge statement, which represent the angles in the Givens rotations,
        in order, that must be applied.
    """

    # QHACK #

    def check(theta1, theta2, theta3, a, b, c, d):

        # 0
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 1
        theta1 = theta1 * -1
        theta2 = theta2
        theta3 = theta3
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 2
        theta1 = theta1
        theta2 = theta2 * -1
        theta3 = theta3
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 3
        theta1 = theta1
        theta2 = theta2
        theta3 = theta3 * -1
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 4
        theta1 = theta1 * -1
        theta2 = theta2 * -1
        theta3 = theta3
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 5
        theta1 = theta1 * -1
        theta2 = theta2
        theta3 = theta3 * -1
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 6
        theta1 = theta1
        theta2 = theta2 * -1
        theta3 = theta3 * -1
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

        # 7
        theta1 = theta1 * -1
        theta2 = theta2 * -1
        theta3 = theta3 * -1
        err1 = np.abs(np.cos(theta1 / 2) * np.cos(theta3 / 2) - a)
        err2 = np.abs(-np.cos(theta1 / 2) * np.sin(theta3 / 2) - d)
        err3 = np.abs(-np.sin(theta1 / 2) * np.cos(theta2 / 2) - b)
        err4 = np.abs(np.sin(theta1 / 2) * np.sin(theta2 / 2) - c)
        if (err1 + err2 + err3 + err4 < 0.0001):
            return theta1, theta2, theta3

    theta1 = 2.0 * np.arcsin(np.sqrt(b * b + c * c))
    theta2 = 2.0 * np.arcsin(c / np.sin(theta1 / 2.0))
    theta3 = 2.0 * np.arcsin(-1.0 * d / np.cos(theta1 / 2.0))
    # check sign
    theta_1, theta_2, theta_3 = check(theta1, theta2, theta3, a, b, c, d)

    return theta_1, theta_2, theta_3
    # QHACK #


if __name__ == "__main__":
    # DO NOT MODIFY anything in this code block
    inputs = sys.stdin.read().split(",")
    theta_1, theta_2, theta_3 = givens_rotations(float(inputs[0]),
                                                 float(inputs[1]),
                                                 float(inputs[2]),
                                                 float(inputs[3]))
    print(*[theta_1, theta_2, theta_3], sep=",")
