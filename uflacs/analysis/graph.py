
from uflacs.analysis.graph_vertices import build_graph_vertices
from uflacs.analysis.graph_symbols import build_graph_symbols
from uflacs.analysis.graph_rebuild import rebuild_expression_from_graph

class Graph2(object):
    def __init__(self):
        self.e2i = {}
        self.V = []
        self.expression_vertices = []
        self.nv = 0
        self.V_shapes = []
        self.V_symbols = None # Crs matrix
        self.total_unique_symbols = 0

def build_graph(expressions, DEBUG=False):

    # Make empty graph
    G = Graph2()

    # Populate with vertices
    G.e2i, G.V, G.expression_vertices = build_graph_vertices(expressions)
    G.nv = len(G.V)

    # Populate with symbols
    G.V_shapes, G.V_symbols, G.total_unique_symbols = \
                build_graph_symbols(G.V, G.e2i, DEBUG)

    if DEBUG:
        assert G.total_unique_symbols == len(set(G.V_symbols.data))

    return G
