import numpy as np
import trimesh
from scipy.sparse import csr_matrix
from collections import *
import heapq
from Flow import *
import queue

ang_eta = 0.3
weight_delta = 0.5
fuzzy_eps = 0.05

class Graph():
    def __init__(self, vertices, g):
        self.V = vertices
        self.graph = self.init_graph(g)

    def init_graph(self, g):
        graph = {}
        for i in range(0, self.V):
            graph[i] = []
            for j in range(g.indptr[i], g.indptr[i+1]):
                graph[i].append((g.indices[j], g.data[j]))
        return graph

    # def min_distance(self, dist, spt_set):
    #     min_val = 1e7
    #     min_index = -1
    #     for v in range(self.V):
    #         if dist[v] < min_val and spt_set[v] == False:
    #             min_val = dist[v]
    #             min_index = v
    #     return min_index

    def dijkstra(self, s):
        dis = defaultdict(lambda:float("inf"))
        dis[s] = 0
        q = [(0,s)]
        vis = set()
        while q:
            _, u = heapq.heappop(q)
            if u in vis: continue
            vis.add(u)
            for v,w in self.graph[u]:
                if dis[v] > dis[u] + w:
                    dis[v] = dis[u] + w
                    heapq.heappush(q,(dis[v],v))
        return dis
        # if src == tgt:
        #     return 0
        # elif self.graph[src, tgt] != 0:
        #     return self.graph[src, tgt]
        # dist = [1e7] * self.V
        # dist[src] = 0
        # spt_set = [False] * self.V
        # for _ in range(self.V):
        #     u = self.min_distance(dist, spt_set)
        #     print(u)
        #     if u == -1:
        #         print("???")
        #         break
        #     spt_set[u] = True
        #     if u == tgt:
        #         break
        #     for v in range(self.V):
        #         if (self.graph[u, v] > 0 and spt_set[v] == False and dist[v] > dist[u] + self.graph[u, v]):
        #             dist[v] = dist[u] + self.graph[u, v]
        # return dist[tgt]
    def get_dist_graph(self):
        dist_list = []
        for i in range(self.V):
            dist_list.append(self.dijkstra(i))
        return dist_list

def floyd(num, graph):
    for i in range(num):
        for j in range(num):
            if i != j and graph[i][j] == 0:
                graph[i][j] = np.inf
    print('here')
    for i in range(num):
        print(i)
        for j in range(num):
            for k in range(num):
                graph[j][k] = min(graph[j][k], graph[j][i] + graph[i][k])
    return graph

def getGeod(idx, verts, norms, convs, face_adj, adj_edges, unsh_edges):
    verts_shared = verts[adj_edges[idx]] # 2, face 0, face 1
    verts_unshared = verts[unsh_edges[idx]] # 2, face 0, face 1
    # origin = np.array(verts_shared[1])
    shared_loc = np.array(verts_shared) - np.array(verts_shared[1])
    unshared_loc = np.array(verts_unshared) - np.array(verts_shared[1]) # (a1, b1, c1), (a2, b2, c2)
    axis = np.array(verts_shared[0]) - np.array(verts_shared[1]) # (a, b, c)
    axis   = axis / np.linalg.norm(axis)
    edge_0 = unshared_loc[0] - shared_loc[1] # shared_loc[1]: origin
    # edge_1 = unshared_loc[1] - shared_loc[1]
    convex = convs[idx]
    norm_0, norm_1 = norms[face_adj[idx]]
    center_0 = sum([shared_loc[0], shared_loc[1], unshared_loc[0]]) / 3
    center_1 = sum([shared_loc[0], shared_loc[1], unshared_loc[1]]) / 3
    lhs_plane = 0 if np.dot(np.cross(axis, edge_0), norm_0) > 0 else 1
    angle = np.arccos(np.dot(norm_0, norm_1))
    if convex == 0:
        angle = -angle
    xn, yn, zn = axis[0], axis[1], axis[2]
    x0, y0, z0 = shared_loc[0][0], shared_loc[0][1], shared_loc[0][2]
    M = xn * x0 + yn * y0 + zn * z0
    c = np.cos(angle)
    s = np.sin(angle)
    rot_mat = np.array(
        [
            [xn * xn * (1 - c) + c, xn * yn * (1 - c) - zn * s, xn * zn * (1 - c) + yn * s, (x0 - xn * M) * (1 - c) + (zn * y0 - yn * z0) * s],
            [xn * yn * (1 - c) + zn * s, yn * yn * (1 - c) + c, yn * zn * (1 - c) - xn * s, (y0 - yn * M) * (1 - c) + (xn * z0 - zn * x0) * s],
            [xn * zn * (1 - c) - yn * s, yn * zn * (1 - c) + xn * s, zn * zn * (1 - c) + c, (z0 - zn * M) * (1 - c) + (yn * x0 - xn * y0) * s],
            [0, 0, 0, 1]
        ]
    ) # whether the 4th term is necessary?
    if lhs_plane == 0:
        rotted_center = np.dot(rot_mat, np.array([center_0[0], center_0[1], center_0[2], 1.0]))[:3]
        geod_dist = np.linalg.norm(rotted_center - center_1)
    else:
        rotted_center = np.dot(rot_mat, np.array([center_1[0], center_1[1], center_1[2], 1.0]))[:3]
        geod_dist = np.linalg.norm(rotted_center - center_0)
    return geod_dist

def getAng(convex, normal1, normal2):
    if convex:
        return ang_eta * (1.0 - np.dot(normal1, normal2))
    else:
        return 1.0 - np.dot(normal1, normal2)

def getProb(idx, idx_seed, seed_list, dist):
    if idx == idx_seed:
        return 1
    s = 0
    for seed in seed_list:
        if dist[idx][seed] == 0:
            return 0
        s += 1 / dist[idx][seed]
    return (1 / dist[idx][idx_seed]) / s

def selectFirstSeed(num_faces, dist):
    key_idx = -1
    key_sum = float("inf")
    for i in range(num_faces):
        d = sum(dist[i].values())
        if d < key_sum:
            key_idx = i
            key_sum = d
    return key_idx

def minDistofSeeds(cand_idx, seed_list, dist):
    return min([dist[cand_idx][seed] for seed in seed_list])

def selectSeeds(num_seeds, num_faces, dist):
    first_seed = selectFirstSeed(num_faces, dist)
    seed_set = set([first_seed])
    for _ in range(num_seeds-1):
        maxmin_dist = 0
        maxmin_idx = -1
        for i in range(num_faces):
            min_dist = minDistofSeeds(i, seed_set, dist)
            if min_dist > maxmin_dist:
                maxmin_dist = min_dist
                maxmin_idx = i
        seed_set.add(maxmin_idx)
    return seed_set

def k_way(num_faces, seed_list, dist):
    res = []
    fuzzy_dict = defaultdict(list)
    for i in range(num_faces):
        max_prob, submax_prob = 0, 0
        max_idx, submax_idx = -1, -1
        for j, seed in enumerate(seed_list):
            prob = getProb(i, seed, seed_list, dist)
            if prob > max_prob:
                max_prob = prob
                max_idx = j
            elif prob > submax_idx:
                submax_prob = prob
                submax_idx = j
        if max_prob - submax_prob >= fuzzy_eps:
            res.append(max_idx)
        else:
            fuzzy_dict[(max_idx, submax_idx)].append(i)
            fuzzy_dict[(submax_idx, max_idx)].append(i)
            res.append(-1)
    return res, fuzzy_dict

def getFuzzyGraph(ang_dist_mat, src_s1, tgt_s2, fuzzy_dict, decomp_res, num_faces):
    # adj_graph: csr_matrix, rows and colomns are face ids
    # src_s1, tgt_s2: INT, id of seeds (depart different regions)
    # fuzzy_dict: dict, {(region_A, region_B):[fuzzy_face_ids]}
    # ang_dict: list, ang_dist of adjacency pairs
    # decomp_res: list, classification of faces
    fuzzy_graph = FlowGraph()
    a1_set = set()
    a2_set = set()
    fuzzy_set = set(fuzzy_dict[(src_s1, tgt_s2)])
    for face in fuzzy_set:
        for i in range(ang_dist_mat.indptr[face], ang_dist_mat.indptr[face+1]):
            c = ang_dist_mat.indices[i]
            if decomp_res[c] == src_s1:
                a1_set.add(c)
                fuzzy_graph.addEdge(c, face, 0)
            elif decomp_res[c] == tgt_s2:
                a2_set.add(c)
                fuzzy_graph.addEdge(face, c, 0)
    for fuzzy_face in fuzzy_set:
        for i in range(ang_dist_mat.indptr[fuzzy_face], ang_dist_mat.indptr[fuzzy_face+1]):
            c = ang_dist_mat.indices[i]
            if c in fuzzy_set:
                fuzzy_graph.addEdge(fuzzy_face, c, 0)

    fuzzy_graph.updateCap(ang_dist_mat)
    return fuzzy_graph



def allo_color(seed_list):
    dic = {seed_list[0]: (50,50,50), seed_list[1]:(50,150,150),seed_list[2]:(50,250,250)}
    return dic


