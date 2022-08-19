#!/usr/bin/env python3

from matplotlib import pyplot as plt
import math
import sys

# Stat to compare refers to statname in node_centrality.tab
# ex) NodeId, Closeness, Betweenes, EigenVector ...

if (len(sys.argv) != 4):
    print("argument count is incorrect!")
    print("usage: python3 test_centrality.py <Stat to compare> <SNAP centrality result> <OUR result>")
    print("example: ./testers/test_centrality.py Betweennes sample-data/centrality/p2p-Gnutella04.snap.tab p2p-Gnutella04_Betweennes.txt")
    exit()

stat_to_compare = sys.argv[1]
base_file = sys.argv[2]
test_file = sys.argv[3]

with open(base_file, mode='r') as f_base, open(test_file, mode='r') as f_test, open("test_output", mode='w') as f_out:
    # base: output from SNAP
    base_lis = []
    stat_idx = -1
    for line in f_base:
        tokens = line.split()
        if tokens[0][0] == '#':
            # is this line a list of stat names?
            if "NodeId" in tokens[0]:
                if stat_to_compare not in tokens:
                    print("stat name not found!")
                    exit()
                else:
                    stat_idx = tokens.index(stat_to_compare)
            continue
        
        if stat_idx == -1:
            print("SNAP's centrality file format is wrong!")
            exit()
        vid = int(tokens[0])
        betweeen = float(tokens[stat_idx])
        base_lis.append((betweeen, vid))
    base_lis.sort(reverse=True)
    
    # test: output from OUR CODE
    test_lis = []
    for line in f_test:
        tokens = line.split()
        if tokens[0][0] == '#':
            continue
        
        vid = int(tokens[0])
        betweeen = float(tokens[1])
        test_lis.append( (betweeen, vid) )
    test_lis.sort(reverse=True)

    test_dict = {}
    for i, (between, vid) in enumerate(test_lis):
        test_dict[vid] = i
    
    f_out.write(f'Comparing Stat: {stat_to_compare}\n\n')
    f_out.write(f'{"SNAP-rank" :<14}{"OUR-rank" :<14}{"vid" :<8}{"SNAP-value" :<16}{"OUR-value" :<16}{"R-Err" :<8}\n')
    for i in range (min(len(base_lis), len(test_lis))):
        (val, vid) = base_lis[i]
        t_idx = test_dict[vid]
        t_val = test_lis[t_idx][0]
        err = str(t_val/val - 1) if (val > 0 and t_val > 0) else "A:" + str(abs(t_val - val))
        l = f'{i :<14}{t_idx :<14}{vid :<8}{val :<16}{t_val :<16}{err :<8}\n'
        f_out.write(l)

