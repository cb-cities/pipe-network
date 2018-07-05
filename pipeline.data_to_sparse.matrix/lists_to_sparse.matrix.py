import json
import numpy as np
import scipy.sparse
from scipy.sparse import csr_matrix
import sys

link_file = 'link_list4graph.json'
link_data = json.load(open(link_file), strict=False)

row = []
col = []
data = []
for dictionary in link_data:
    row = np.append(row, dictionary['link_id'])
    row = np.append(row, dictionary['link_id'])
    col = np.append(col, dictionary['start_node'])
    col = np.append(col, dictionary['end_node'])
    if dictionary['start_node'] < dictionary['end_node']:
        data = np.append(data, [-1,1])
    else:
        data = np.append(data, [1,-1])

col = col.astype(int)
row = row.astype(int)
#csr_matrix((data,(row,col)))
m = csr_matrix((data,(row,col)))

#np.save('sparse_matrix.npy', m)
scipy.sparse.save_npz('sparse_matrix.npz', m)
#A = scipy.load('sparse_matrix.npy')[()]