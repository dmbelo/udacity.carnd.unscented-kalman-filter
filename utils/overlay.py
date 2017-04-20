"""
Overlay plot
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt

def overlay(file_name):

    column_names = ['px_est', 'py_est', 'px_meas', 'py_meas', 'px_gt', 'py_gt'] 

    with open(file_name) as f:
        df = pd.read_table(f, index_col=False, sep='\t', header=None,
                        names=column_names, lineterminator='\n');

    plt.plot(df.px_gt, df.py_gt, 'o', df.px_est, df.py_est, '+')
    plt.legend(['Ground Truth', 'Estimation'])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

if __name__ == "__main__":

    file_name = sys.argv[1]
    overlay(file_name)