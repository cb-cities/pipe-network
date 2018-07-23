import json

pipeline_file='wPressurizedMain.json'
pipe_data = json.load(open(pipeline_file), strict=False)

p_f = pipe_data['features']

node_list = []
for feature in p_f:
    coord1 = feature['geometry']['coordinates'][0]
    node_list.append((round(coord1[0],12), round(coord1[1],12)))
    coord2 = feature['geometry']['coordinates'][-1]
    node_list.append((round(coord2[0],12), round(coord2[1],12)))

#    for coord in feature['geometry']['coordinates']:
#        node_list.append((round(coord[0],12), round(coord[1],12)))

node_set = set(node_list)


node_list2 = list(node_set)

node_id = 0
node_dict = {}
for node in node_list2:
    node_dict[node]=node_id
    node_id += 1

node_dict2 = {str(k):v for k, v in node_dict.items()}


with open('node_temp_bz247.json', 'w') as outfile:
    json.dump(node_dict2, outfile, indent=2)

node_list4graph = []

for key,value in node_dict.items():
    temp_node_dict = {'node_id': value, 'node_coord_x': key[0], 'node_coord_y': key[1]}
    node_list4graph.append(temp_node_dict)

link_list4graph=[]


link_id = 0
for link in p_f:
    temp_start_x = round(link['geometry']['coordinates'][0][0],12)
    temp_start_y = round(link['geometry']['coordinates'][0][1],12)
    temp_end_x = round(link['geometry']['coordinates'][-1][0],12)
    temp_end_y = round(link['geometry']['coordinates'][-1][1],12)
    temp_link_dict = {'link_id': link_id, 'start_node': node_dict[(temp_start_x,temp_start_y)], 'end_node': node_dict[(temp_end_x, temp_end_y)]}
    link_id += 1
    link_list4graph.append(temp_link_dict)

with open('node_list4graph_bz247.json', 'w') as outfile:
    json.dump(node_list4graph, outfile, indent=2)

with open('link_list4graph_bz247.json', 'w') as outfile:
    json.dump(link_list4graph, outfile, indent=2)
