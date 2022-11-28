import numpy as np
import trimesh
from scipy.sparse import csr_matrix
from collections import *
from Flow import *
import queue

# mesh = trimesh.load_mesh("test.ply")
# verts = mesh.vertices
# norms = mesh.face_normals
# faces = mesh.faces
# face_adjacency = mesh.face_adjacency
# face_adjacency_convex = mesh.face_adjacency_convex
# face_adjacency_edges = mesh.face_adjacency_edges
# face_adjacency_unshared = mesh.face_adjacency_unshared
# print(verts)
# print(faces)
# print(face_adjacency)
# print(face_adjacency_unshared)

# a = np.random.randint(0,10, size=(3,3))
# b = np.array([1, 2, 3])
# print(a)
# print(b)
# print(np.dot(a, b))

# row = np.array([0, 1, 2, 0])
# col = np.array([0, 1, 1, 0])
# data = np.array([1, 2, 4, 8])
# csr = csr_matrix((data, (row, col)), shape=(3, 3))
# print(csr + csr.getH())

# dist = [[0 for _ in range(3)] for __ in range(3)]
# print(dist)

# a = [5,3,2,3,1]
# print(min(a))

# mesh = trimesh.load_mesh("bunny.ply")
# v = mesh.vertices
# f = mesh.faces

# v = np.array(v)
# f = np.array(f)

# mesh = trimesh.Trimesh(vertices = v, faces = f)

# mesh = trimesh.Trimesh(vertices = v, faces = f, process = False)

# for facet in range(2000):
#     mesh.visual.face_colors[facet] = trimesh.visual.random_color()
# mesh.show()

# print(mesh.faces[0])

# fuzzy_graph = defaultdict(list)
# print(fuzzy_graph[3])

# g = FlowGraph()
# g.addEdge(0, 4, 1000)
# g.addEdge(4, 3, 20)
# g.addEdge(4, 2, 30)
# g.addEdge(2, 1, 10)
# g.addEdge(2, 3, 10)
# g.addEdge(1, 3, 30)
# g.addEdge(3, 4, 20)
# g.addEdge(2, 4, 30)
# g.addEdge(1, 2, 10)
# g.addEdge(3, 2, 10)
# g.addEdge(3, 1, 30)
# g.addEdge(3, 5, 1000)
# print(g.maxFlow(0, 5))
# print(g.BFS(4))

# q = queue.Queue()
# q.put(1)
# q.put(3)
# q.put(2)
# print(q.qsize())
# print(q.get())
# print(q.qsize())
# print(q.get())
# print(q.qsize())


# a = set()
# a.add(0)
# a.add(2)
# a.add(1)
# print(len(a))
# for i in a:
#     print(i)
# print(3 in a)
# print(3 not in a)