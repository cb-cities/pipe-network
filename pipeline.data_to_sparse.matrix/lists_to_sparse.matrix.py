import json
import numpy as np
from scipy.sparse import csr_matrix

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

csr_matrix((data,(row,col)))
 
m = csr_matrix((data,(row,col)))
np.save('sparse_matrix.npy', m)