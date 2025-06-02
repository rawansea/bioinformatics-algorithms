# Rawan Abu Alkhayr
 

# Sources: 
# https://www.youtube.com/watch?v=f5kgmqcwb8M
# https://github.com/TatyanaV/Genome_Sequencing_Bioinformatics_II/blob/master/8.StringReconstructionProblem.py
# https://github.com/quevedito2/Bioinformatics-coursera/blob/master/22_StringReconstruction.py


from collections import defaultdict

# function to build de Bruijn graph
def build_de_bruijn_graph(reads, k):
    edges = []
    nodes = set()
    for read in reads:
        for i in range(len(read) - k + 1):
            prefix = read[i:i + k - 1]
            suffix = read[i + 1:i + k]
            edges.append((prefix, suffix))
            nodes.add(prefix)
            nodes.add(suffix)
    return nodes, edges

class DeBruijnGraph:
    def __init__(self, reads, k):
        self.graph = defaultdict(list)
        self.in_deg = defaultdict(int)
        self.out_deg = defaultdict(int)
        _, edges = build_de_bruijn_graph(reads, k)
        for u, v in edges:
            self.graph[u].append(v)
            self.out_deg[u] += 1   # Track out-degrees
            self.in_deg[v] += 1    # Track in-degrees
            
    # Perform Eulerian walk using Hierholzer’s algorithm
    def find_eulerian_path(self):
        # Find starting node
        start = None
        for node in self.graph:
            if self.out_deg[node] > self.in_deg[node]:
                start = node
                break
        if start is None:
            start = next(iter(self.graph))
        
        stack = [start]
        path = []
        while stack:
            current = stack[-1]
            if self.graph[current]:
                next_node = self.graph[current].pop()  # Explore edge
                stack.append(next_node)
            else:
                path.append(stack.pop())
        return path[::-1]
        
# Reconstruct the original string
def reconstruct_superstring(path):
    if not path:
        return ""
    return path[0] + ''.join(node[-1] for node in path[1:])


# ----------------------------

# User input
k = int(input("Enter the value of k: "))
reads = input("Enter all reads separated by spaces: ").strip().split()

# Build graph and reconstruct genome
G = DeBruijnGraph(reads, k)
path = G.find_eulerian_path()
superstring = reconstruct_superstring(path)

print("Assembled genome:", superstring)
