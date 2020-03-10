from scipy.spatial import Delaunay
import numpy as np


class Shape:
    def __init__(self, coords):
        self.coords = coords

    def alphashape(self, alpha, computer_inner=False):
        if len(self.coords) < 4:
            raise Exception("At least four points are needed for alpha shape computation")

        tri = Delaunay(self.coords)
        edges = set()

        def add_edge(x, y):
            if (x, y) in edges or (y, x) in edges:
                if not computer_inner:
                    edges.remove((y, x))
            else:
                edges.add((x, y))

        def calc_circumradius(x, y, z):
            px = self.coords[x]
            py = self.coords[y]
            pz = self.coords[z]

            a = np.sqrt((px[0] - py[0]) ** 2 + (px[1] - py[1]) ** 2)
            b = np.sqrt((py[0] - pz[0]) ** 2 + (py[1] - pz[1]) ** 2)
            c = np.sqrt((pz[0] - px[0]) ** 2 + (pz[1] - px[1]) ** 2)

            s = (a + b + c) / 2.0

            area = np.sqrt(s * (s - a) * (s - b) * (s - c))
            return a * b * c / (4.0 * area)

        for ix, iy, iz in tri.vertices:
            circumradius = calc_circumradius(ix, iy, iz)
            if circumradius < alpha:
                add_edge(ix, iy)
                add_edge(iy, iz)
                add_edge(iz, ix)

        return edges
