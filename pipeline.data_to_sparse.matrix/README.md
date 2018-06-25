# pipeline.data_to_sparse.matrix

These are the python script which converts the pipeline data to sparse matrix which represents the connection between nodes and links. 

## pipeline.data_to_lists.py
This is the script for preparing the files to create a sparse matrix.

* Email me and get the original data file, "wPressurizedMain.json" which is written in [pipeline.data_to_lists.py](pipeline.data_to_lists.py).

* Executing the script, you will get these json files;
  - 'node_temp.json': a dictionary which has coordinates of nodes as keys and node ID as values.
  - 'node_list4graph.json': a list which has dictionaries in it. Each dictionary has a information about each node, which is node ID, X(longtitude) and Y(latitude) coordinates of the node as keys. 
  - 'link_list4graph.json': a list which has dictionaries in it. Each dictionary has a information about each link, which is link ID and the start and end node ID of the link as keys.

## lists_to_sparse.matrix.py
This is the script for creating a sparse matrix.

