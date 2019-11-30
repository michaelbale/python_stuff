# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 01:21:34 2018

@author: Michael
"""

import numpy as np
import sys

def get_oci(sorted_data):
    num_sites = np.sum(sorted_data)
    num_uis = len(sorted_data)
    rel_abun = sorted_data / num_sites
    cum_rel_abun_list = list()
    cum_rel_abun_list.append(rel_abun[0])
    for x in range(1, num_uis):
        cum_rel_abun_list.append(cum_rel_abun_list[x-1] + rel_abun[x])
    cum_rel_abun = np.array(cum_rel_abun_list)
    oligoclonality_index = 2 * (np.sum(cum_rel_abun/num_uis) - 0.5)
    return oligoclonality_index
    


def main(argv):
    with open(argv[0], "r") as fh:
        content = fh.readlines()
    content = [x.strip() for x in content]
    content = [int(x) for x in content]
    data = np.array([x for x in content if x > 0])
    sorted_data = -np.sort(-data)
    print(get_oci(sorted_data))



if __name__ == '__main__':
    main(sys.argv[1:])