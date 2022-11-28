import numpy as np
import queue
from collections import *

class FlowEdge():
    def __init__(self, fr, to, cap, flow):
        self.fr = fr
        self.to = to
        self.cap = cap
        self.flow = flow

class FlowGraph():
    def __init__(self):
        # self.V = 0
        self.edges = []
        self.graph = defaultdict(list)
        # self.a = defaultdict(int)
        self.p = {}
    def addEdge(self, fr, to, cap):
        self.edges.append(FlowEdge(fr, to, cap, 0))
        self.edges.append(FlowEdge(to, fr, 0, 0))
        m = len(self.edges)
        self.graph[fr].append(m - 2)
        self.graph[to].append(m - 1)
    def maxFlow(self, src, tgt):
        # print(self.graph)
        flow = 0
        # v = len(self.graph)
        while True:
            a = defaultdict(int)
            Q = queue.Queue()
            Q.put(src)
            a[src] = float("inf")
            while Q.empty() == False:
                x = Q.get()
                for idx in self.graph[x]:
                    e = self.edges[idx]
                    if a[e.to] == 0 and e.cap > e.flow:
                        self.p[e.to] = idx
                        a[e.to] = min(a[x], e.cap - e.flow)
                        Q.put(e.to)
                        # print(e.to, a[e.to])
                # print(a[1], a[2], a[3], a[4])
                if a[tgt] != 0:
                    break
            if a[tgt] == 0:
                break
            u = tgt
            while True:
                if u == src:
                    break
                self.edges[self.p[u]].flow += a[tgt]
                self.edges[self.p[u] ^ 1].flow -= a[tgt]
                u = self.edges[self.p[u]].fr
            # print(a[tgt])
            flow += a[tgt]
        return flow
    def BFS(self, src):
        searched_v = set()
        Q = queue.Queue()
        Q.put(src)
        searched_v.add(src)
        while Q.empty() == False:
            x = Q.get()
            for idx in self.graph[x]:
                e = self.edges[idx]
                if idx % 2 == 0 and e.cap > e.flow and e.to not in searched_v:
                    searched_v.add(e.to)
                    Q.put(e.to)
        not_searched = []
        for v in self.graph:
            if v not in searched_v:
                not_searched.append(v)
        return list(searched_v), not_searched
    def updateCap(self, ang_mat):
        # Not including src and tgt
        ang_list = [ang_mat[e.fr, e.to] for e in self.edges]
        ang_list = np.array(ang_list)
        ang_mean = ang_list.mean()
        for i in range(len(self.edges)):
            if i % 2 == 0:
                self.edges[i].cap = 1 / (1 + ang_mat[self.edges[i].fr, self.edges[i].to] / ang_mean)

