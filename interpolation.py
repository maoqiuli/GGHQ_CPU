import string
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


def Interpolation(file, tar_QPS):

    with open(file, "r") as f:
        data = f.readlines()


    data = [data[i].split("\t") for i in range(len(data))]

    recall_ind = [i for i in range(len(data[0])) if data[0][i] == "R@10"][0]
    QPS_ind = [i for i in range(len(data[0])) if data[0][i] == "QPS"][0]

    recall = [float(data[i][recall_ind]) for i in range(1, len(data))]
    QPS = [float(data[i][QPS_ind]) for i in range(1, len(data))]

    tar_recall = np.interp(tar_QPS, QPS[::-1], recall[::-1])

    return tar_recall