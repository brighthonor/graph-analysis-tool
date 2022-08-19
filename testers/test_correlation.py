#!/usr/bin/env python3

import sys

# file[0]: file name, file[1]: column
def get_data(file):
    fname, column = file[0], int(file[1])
    res = {}
    with open(fname, mode = 'r') as f:
        for line in f:
            cols = line.split()
            if not cols[0].isnumeric():
                continue # numeric id only
            res[int(cols[0])] = float(cols[column])
    return res

# data[0], data[1]: dictionary with same keys
# returns tuple (a, b, r^2)
def calc_correlation(data):
    sumx, sumx2, sumy, sumy2, sumxy, n = 0, 0, 0, 0, 0, 0
    for key in data[0]:
        assert key in data[1], f'Error: key {key} does not exist in file 2'
        x, y = data[0][key], data[1][key]
        sumx += x
        sumy += y
        sumx2 += x * x
        sumy2 += y * y
        sumxy += x * y
        n += 1

    b = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx)
    a = (sumy - b * sumx) / n
    r2 = b * (n * sumxy - sumx * sumy) / (n * sumy2 - sumy * sumy)
    return a, b, r2

# input: list of 2-tuples sorted by the second item of each tuple
# output: dictionary that contains the rank of each item
# example:
#     input: [(1, 3), (5, 4), (3, 4), (0, 6)]
#     output: {1: 0, 5: 1.5, 3: 1.5, 0: 3}
def key_to_rank(data):
    ranks = dict() # calculate rank ranges (to handle ties correctly)
    for i in range(len(data)):
        if data[i][1] in ranks:
            ranks[data[i][1]][1] = i
        else:
            ranks[data[i][1]] = [i, i]
    return {k[0]: (ranks[k[1]][0] + ranks[k[1]][1]) * 0.5 for k in data}

# Spearman's rank correlation coefficient (squared)
# returns single value between 0 ~ 1
def calc_spearman_rank_correlation(data):
    # convert dictionaries to lists of 2-tuples (key, value) sorted by value
    sorted_items = [sorted(a.items(), key = lambda x: x[1], reverse = True) for a in data]
    # convert values to ranks
    item_ranks = [key_to_rank(x) for x in sorted_items]
    return calc_correlation(item_ranks)[2]

if len(sys.argv) != 5:
    print('usage: test_correlation.py <file 1> <column> <file 2> <column>')
    print('example:')
    print('$ ./main -uf sample-data/p2p-Gnutella08.snap')
    print('$ ./testers/test_correlation.py sample-data/centrality/p2p-Gnutella08.snap.tab 4')
    print('      p2p-Gnutella08_EigenCentrality.txt 1')
    exit(1)

files = [sys.argv[1:3], sys.argv[3:5]]

data = [get_data(f) for f in files]

a, b, r2 = calc_correlation(data)

rho2 = calc_spearman_rank_correlation(data)

print('y = a + b x')
print(f'a = {a}')
print(f'b = {b}')
print(f'r^2 = {r2}')
print(f"Spearman's rank correlation coefficient:\nrho^2 = {rho2}")
