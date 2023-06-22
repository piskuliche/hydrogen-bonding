#!/usr/bin/env python3

import numpy as np
import math

def read_frame(f, don_atoms, acc_atoms):
    line = f.readline().strip().split()
    donors = {0:[],1:[]}
    acceptors = {0:[],1:[]}
    hbonds = []
    for i in [0,1]:
        line = f.readline().strip()
        num_hbonds = int(line)
        hbonds.append(num_hbonds)
        for j in range(num_hbonds):
            line = f.readline().strip().split()
            donors[i].append(math.ceil(int(line[0])/don_atoms))
            acceptors[i].append(math.ceil(int(line[1])/acc_atoms[i]))
    return hbonds

if __name__ == "__main__":

    with open("hbonding.out",'r') as f:
        hbonds = {0:[],1:[]}
        for i in range(500):
            print(i)
            hbtmp = read_frame(f, 2, [9,2])
            hbonds[0].append(hbtmp[0])
            hbonds[1].append(hbtmp[1])
        print(np.average(hbonds[0]))
        print(np.average(hbonds[1]))