from tools import *

mesh = trimesh.load_mesh("bunny.ply")

face_list = mesh.faces
num_faces = face_list.shape[0]
face_normals = mesh.face_normals
vertices = mesh.vertices

face_adjacency = mesh.face_adjacency
face_adjacency_convex = mesh.face_adjacency_convex
face_adjacency_edges = mesh.face_adjacency_edges
face_adjacency_unshared = mesh.face_adjacency_unshared

geod_dists = []
ang_dists = []

for idx, (face_idx_1, face_idx_2) in enumerate(face_adjacency):
    angular_dist = getAng(face_adjacency_convex[idx], face_normals[face_idx_1], face_normals[face_idx_2])
    ang_dists.append(angular_dist)
    geodesic_dist = getGeod(idx, vertices, face_normals, face_adjacency_convex, face_adjacency, face_adjacency_edges, face_adjacency_unshared)
    geod_dists.append(geodesic_dist)
ang_dists = np.array(ang_dists)
geod_dists = np.array(geod_dists)

rows = [pair[0] for pair in face_adjacency]
cols = [pair[1] for pair in face_adjacency]
weight = weight_delta * geod_dists / geod_dists.mean() + (1.0 - weight_delta) * ang_dists / ang_dists.mean()

weight_mat = csr_matrix((weight, (rows, cols)), shape=(num_faces, num_faces))
ang_dist_mat = csr_matrix((ang_dists, (rows, cols)), shape=(num_faces, num_faces))

weight_mat = weight_mat + weight_mat.getH()

graph = Graph(num_faces, weight_mat)
dist = graph.get_dist_graph()

seed_set = selectSeeds(4, num_faces, dist)
seed_list = list(seed_set)

res = k_way(num_faces, seed_list, dist)

color_dic = {seed_list[0]: (240,100,160,150), seed_list[1]:(50,150,150,150),seed_list[2]:(250,250,100,150),seed_list[3]:(150,100,250,150)}

for face in range(num_faces):
    mesh.visual.face_colors[face] = color_dic[seed_list[res[face]]]
mesh.show()
