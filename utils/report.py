"""
Overlay plot
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def overlay(file_name):

    column_names = ['time', 'px_est', 'py_est', 'v_est', 'a_yaw_est', 'n_yaw_est', 'px_meas', 'py_meas', 'px_gt', 'py_gt', 'vx_gt', 'vy_gt', 'NIS'] 

    with open(file_name) as f:
        df = pd.read_table(f, index_col=False, sep='\t', header=None,
                           names=column_names, lineterminator='\n');

    vx_est = df.v_est * np.cos(df.a_yaw_est)
    vy_est = df.v_est * np.sin(df.a_yaw_est)

    plt.figure()
    plt.plot(df.px_gt, df.py_gt, df.px_meas, df.py_meas, 'o', df.px_est, df.py_est, 'k+')
    plt.legend(['Ground Truth', 'Measurement', 'Estimation'])
    plt.xlabel('X Position')
    plt.ylabel('Y Position')

    plt.figure()
    plt.plot(df.vx_gt, df.vy_gt, vx_est, vy_est, 'k+')
    plt.legend(['Ground Truth', 'Estimation'])
    plt.xlabel('X Velocity')
    plt.ylabel('Y Velocity')

    x = df.time.iloc[[0, -1]]
    y = 7.81 * np.ones([2, 1])
    plt.figure()
    plt.plot(df.time, df.NIS, x, y)
    plt.ylabel('NIS')

    plt.figure()
    plt.subplot(221)
    plt.plot(df.time, df.px_gt, df.time, df.px_est, '+')
    plt.xlabel('Time')
    plt.ylabel('Px')

    plt.subplot(222)
    plt.plot(df.time, df.py_gt, df.time, df.py_est, '+')
    plt.xlabel('Time')
    plt.ylabel('Py')

    plt.subplot(223)
    plt.plot(df.time, df.vx_gt, df.time, vx_est, '+')
    plt.xlabel('Time')
    plt.ylabel('Vx')
    
    plt.subplot(224)
    plt.plot(df.time, df.vy_gt, df.time, vy_est, '+')
    plt.xlabel('Time')
    plt.ylabel('Vy')
    # plt.subplot(121)
    # plt.plot(df.time, df.px_gt, df.time, df.px_est)
    # plt.subplot(122)
    # plt.plot(df.time, error**2)

    plt.show()

if __name__ == "__main__":

    file_name = sys.argv[1]
    overlay(file_name)