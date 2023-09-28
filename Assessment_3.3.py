import heapq
import math
import pandas as pd

class Trie:
    def __init__(self, char, freq, left=None, right=None):
        self.char = char
        self.freq = freq
        self.left = left
        self.right = right
    
    def __lt__(self, other):
        return self.freq < other.freq
    
    

class Coordinates:
    def __init__(self, x, y):
            self._x = x
            self._y = y
            

    def __str__(self):
        return f"Coordinates: ({self._x}, {self._y})"
            
# define the water network class
class WaterNetwork:    
    def __init__(self, size=0):
        self._size = size
        self._node_list = {}
        self._edge_list = {}
        
        
    
    def __str__(self):
        return f"Number of nodes: {self._size}"
    
    # define the node class inside the network
    class Node:
        def __init__(self, key=None, x = 0, y = 0):
            self._key = key
            self._coordinates = Coordinates(x, y)
            self._type = []
            self._link = []
            
            
        def __str__(self):
            return f"Key: {self._key}, {self._coordinates}, Type: {self._type}, Link: {self._link}"

        # Calculate distance between 2 nodes
        def calculate_distance(self, dest):
            return math.sqrt((self._coordinates._x - dest._coordinates._x) ** 2 + (self._coordinates._y - dest._coordinates._y) ** 2)
        
        # Find the neighbours of a vertex
        def neighbour_of_vertex(self):
            return self._link


    def __len__(self):
        return self._size
    
    # There are 2 options (1: 1-way (directed) edge list, 2: 2-way (undirected) edge list)
    def create_edge_list(self, option):
        # Initialize pairs of which keys ranging from 1 to the number of the nodes and values are empty list
        self._edge_list = {i:[] for i in range(1, self._size + 1)}
        # Loop through the edge list to assign the them into the self._edge_list dictionary
        for i in range(1, self._size + 1):  
            for j in range(len(self._node_list[i]._link)):
                src_key = i
                dest_key = self._node_list[src_key]._link[j]
                src_node = self._node_list[src_key]
                try:
                    dest_node = self._node_list[dest_key]
                except:
                    continue
                dist = src_node.calculate_distance(dest_node)
                    # This dictionary structure: Source:[Destination, Cost]
                self._edge_list[src_key].append([dest_key, dist])     
                # if the edge is 2-way
                if option == 2:
                    self._edge_list[dest_key].append([src_key, dist])
        return self._edge_list
    
    
    # This function uses dfs traverse algorithm to detect cycles
    def dfs_detect_cycles(self, vertex, visited, pre_vertex, path, cycles):
        # Labels every node as not visited
        visited[vertex] = True
        # Start with selected vertex
        path.append(vertex)
        # Loop through all its neighbour
        for neighbour in self._node_list[vertex]._link:
            # Ignore the vertex 0
            if neighbour == 0:
                continue
            # Visit this neighbour if it has not been visited yet
            if not visited[neighbour]:
                self.dfs_detect_cycles(neighbour, visited, vertex, path, cycles)
            # If encountering a met neighbour
            if neighbour != pre_vertex and neighbour in path:
                # Get the index of this neighbour in the path
                cycle_start = path.index(neighbour)
                # The cycle starts from this neighbour toward the end of the path
                cycle = path[cycle_start:]
                # Append the new cycle to cycle list
                cycles.append(cycle)
                # Print the new cycle
                print("Cycle:", " -> ".join(map(str, cycle + [cycle[0]])))
                
        # After visiting all neighbours, switch to a different node to explore
        path.pop()
    
    
    # This function helps to create a dict of edges between headwaters of this graph (assume it is a complete graph)
    def create_edge_complete_graph(self):
        # Create a dictionary of edges
        edge_list = {i:[] for i in range(self._size + 1)}
        # Loop through all nodes in node list
        for src in self._node_list:
            for dest in self._node_list:
                if src == dest:
                    continue
                else:
                    # Append a headwater, its neighbour, and the distance to the dictionary
                    if 'headwater' in self._node_list[src]._type and 'headwater' in self._node_list[dest]._type:
                        dist = self._node_list[src].calculate_distance(self._node_list[dest])
                        edge_list[src].append([dest, dist])
                        edge_list[dest].append([src, dist])
        return edge_list
                    
    
    # This function helps to find headwaters in a specific area
    def headwaters_in_range(self, coord1, coord2):
        headwaters = {}
        # iterate through each node and check if that node in the specific range and its type is headwater then append to the headwater dictionary specified above
        for node in self._node_list:
            type = self._node_list[node]._type 
            x = self._node_list[node]._coordinates._x 
            y = self._node_list[node]._coordinates._y 
            if 'headwater' in type and x > coord1._x and x < coord2._x and y > coord1._y and y < coord2._y:
                headwaters[node] = self._node_list[node]
                
        return headwaters   # return all the headwaters met the conditions
    

    # This function helps create a Prim minimum spanning tree from the adjacency list created by the function "weighted_edge_list()"
    def prim(self, adj, top_right_headwater):
        # Initialize an empty dictionary to store the minimum spanning tree
        prim_mst = []
        # Initialize a set to keep track of visited vertices
        visit = set()
        # Get the keys (vertices) from the adjacency list
        # Choose the first key as the starting point
        first_key = top_right_headwater
        # Create a min-heap to prioritize edges based on their min weights
        minH = [[0, first_key]]
        # Continue until all vertices have been visited
        while len(visit) < len(adj):
            # Extract the vertex with the minimum edge cost from the min-heap
            cost, i = heapq.heappop(minH)
            # If the vertex is already visited, skip it
            if i in visit:
                continue
            # Add the vertex to the minimum spanning tree
            prim_mst.append(i)
            # Mark the vertex as visited
            visit.add(i)
            # Iterate over the neighboring vertices and their edge costs
            for neighbour, neighbour_cost in adj[i]:
                # If the neighboring vertex has not been visited, add it to the min-heap
                if neighbour not in visit:
                    # print("i", i)
                    # print(f"Neighbour: {neighbour} and neighbour cost: {neighbour_cost}")
                    heapq.heappush(minH, [neighbour_cost, neighbour])
        # Return the Prim minimum spanning tree
        return prim_mst

            
    def edge_list_in_range(self, coord1, coord2, edge_list):   
        # Get headwaters in range  
        headwaters_in_range = self.headwaters_in_range(coord1, coord2)
        edge_list_in_range = {i:[] for i in range(self._size + 1)}
        for src in edge_list:
            if src in headwaters_in_range:
                for dest in edge_list[src]:
                    if dest[0] in headwaters_in_range:
                        edge_list_in_range[src].append(dest)
                        edge_list_in_range[dest[0]].append([src, dest[1]])
        
        filtered_dict = {key: value for key, value in edge_list_in_range.items() if len(value) != 0}
        return filtered_dict
    
    
    # This function helps find the top right headwater 
    def top_right_headwater(self, coord1, coord2, adj):
        min_dist_to_top_right = float('inf')
        top_right_node = None
        
        for headwater in adj:
            # Calculate distance from the headwater to the top right point
            dist_to_top_right = self._node_list[headwater].calculate_distance(self.Node(x= coord2._x, y=coord1._y))
            # Finding the closest headwater to the top right point
            if dist_to_top_right < min_dist_to_top_right:
                top_right_node = headwater
                min_dist_to_top_right = dist_to_top_right
        return top_right_node           
    
    # This function will return the prim mst  
    def shortest_flight_path_search(self, coord1, coord2):
        adj = self.create_edge_complete_graph()
        adj_in_range = self.edge_list_in_range(coord1, coord2, adj) 
        top_right_headwater = self.top_right_headwater(coord1, coord2, adj_in_range)
        prim = self.prim(adj_in_range, top_right_headwater)
        return prim
    
        
    # This function helps find the original sources
    def bfs_find_sources(self, initial_vertex):
        # Initialize a list to keep track of visited vertices
        visited = [False] * self._size
        
        # Create an edge list based on a certain option (option=2)
        self.create_edge_list(option=2)
        
        # Initialize variables to track source discovery
        source_found = False
        sources = []
        
        # Initialize a queue for the Breadth-First Search (BFS)
        queue = []
        queue.append(initial_vertex)
        
        # Start the breadth-first search to search for the nearest headwater
        while(len(queue) > 0):
            # Get the current node from the queue
            curr_vertex = queue.pop()
            
            # Find its neighbors by iterating through edges
            for edge in self._edge_list[curr_vertex]:
                neighbour = edge[0]
                
                # If the neighbor has not been visited, add it to the queue
                if not visited[neighbour]:
                    queue.append(neighbour)
                
                # Check if the neighbor is a "headwater" source
                if "headwater" in self._node_list[neighbour]._type:
                    source_found = True
                    sources.append(neighbour)
            
            # If a source is found, return the list of sources
            if source_found:
                return sources
                
            
    # def dfs_detect_sources(self, vertex, visited, concentration, max_weight, weight, sources:[], pre_vertex, edge_list):
    #     visited[vertex] = True
    #     if vertex == 56 or vertex == 57:
    #         pass
    #     if "junction" in self._node_list[vertex]._type:
    #         weight += concentration[vertex]
    #     elif "headwater" in self._node_list[vertex]._type and weight == max_weight:
    #         sources.append(vertex)
    #         weight -= concentration[pre_vertex]
            
    #     # for neighbour in self._node_list[vertex]._link:
    #     for neighbour in edge_list[vertex]:
    #         v = neighbour[0]
    #         if v == 0:
    #             continue
    #         if not visited[v]:
    #             self.dfs_detect_sources(v, visited, concentration, max_weight, weight, sources, vertex, edge_list)
         
        
    # def chemical_sources(self, sequence_of_junctions:tuple):
    #     concentration = {i : 0 for i in range(self._size)}
    #     max_weight = 0
    #     min_conc = float("inf")
    #     min_junc = None
    #     for junction in sequence_of_junctions:
    #         concentration[junction[0]] = junction[1]
    #         max_weight += junction[1]
    #         if junction[1] < min_conc:
    #             min_junc = junction[0]
    #             min_conc = junction[1]
    #     # print(concentration)
    #     edge_list = self.create_edge_list(option=2)
    #     visited = [False] * (self._size + 1)
    #     sources = []
    #     self.dfs_detect_sources(min_junc, visited, concentration, max_weight, 0, sources, None, edge_list)
             
    #     return sources
     
    # This function helps find the largest-concentration junction and return sources
    def chemical_sources(self, sequence_of_junctions:tuple):
        max_conc = sequence_of_junctions[0][1]
        max_junc = sequence_of_junctions[0][0]
        for junction in sequence_of_junctions:
            if junction[1] > max_conc:
                max_conc = junction[1]
                max_junc = junction[0]
                
        return self.bfs_find_sources(max_junc)
# This function helps read location data from the csv file
def read_csv_file(file_path):
    location_list = pd.read_csv(file_path)
    file_size = len(location_list)
    first_row = location_list.iloc[0]
    water_network = WaterNetwork()
    base_key = first_row['Node']
    
    for index in range(0, file_size): 
        row = location_list.iloc[index]
        key = row['Node']
        x = row['x']
        y = row['y']
        type = row['type']
        link = row['linked']
        if key == 1:
            water_network._node_list[key] = water_network.Node(key, x, y)
        if key != base_key:
            base_key = key
            water_network._node_list[key] = water_network.Node(key, x, y)
        water_network._node_list[key]._type.append(type)
        water_network._node_list[key]._link.append(link)
        water_network._size = len(water_network._node_list)
        # Create 1-way edge list
        water_network.create_edge_list(option=1)
    return water_network

# This function helps red the names of rivers and creeks
def read_river_names_file(file_path):
    name_list = pd.read_csv(file_path)
    file_size = len(name_list)
    name_text = ''
    for i in range(file_size):
        row = name_list.iloc[i][0]
        name_text = name_text + row
    return name_text


def get_encoding_trie(text):
    if len(text) == 0:
        return
    # Find the frequency of the character in the text
    freq = {ch: text.count(ch) for ch in set(text)}
    # Create a list of trie node (key, value)
    pq = [Trie(key, value) for key, value in freq.items()] 
    # Convert the list into a heap 
    heapq.heapify(pq)
    while len(pq) > 1:
        # Pop the two nodes with the lowest frequencies from the heap
        left, right = heapq.heappop(pq), heapq.heappop(pq)
        # Create a new node with a frequency equal to the sum of the frequencies of the two children.
        new_freq = left.freq + right.freq
        # Push the new node back into the heap.
        heapq.heappush(pq, Trie(None, new_freq, left, right))
    root = pq[0]
    return root


# This function generates binary codes for characters in the encoding trie.
def generate_binary_codes(root):
    # Internal recursive function to traverse the encoding trie.
    def _generate_binary_codes(node, current_code, result):
        # If the current node represents a character, add its binary code to the result.
        if node.char is not None:
            result[node.char] = current_code
            return
        
        # Recursively traverse the left subtree, appending '0' to the current binary code.
        if node.left:
            _generate_binary_codes(node.left, current_code + '0', result)
        
        # Recursively traverse the right subtree, appending '1' to the current binary code.
        if node.right:
            _generate_binary_codes(node.right, current_code + '1', result)

    # Initialize an empty dictionary to store binary codes.
    binary_codes = {}
    
    # Start the traversal from the root of the encoding trie with an empty current_code.
    _generate_binary_codes(root, '', binary_codes)
    
    return binary_codes  # Return the dictionary of binary codes.

# This function encodes an input string using binary codes from an encoding trie.
def code_string(root, input_string):
    # Generate binary codes for characters in the encoding trie.
    binary_codes = generate_binary_codes(root)
    
    # Initialize an empty string to store the encoded output.
    encoded_string = ""

    # Iterate through each character in the input string.
    for char in input_string:
        # Check if the character has a corresponding binary code.
        if char in binary_codes:
            # Append the binary code of the character to the encoded string.
            encoded_string += binary_codes[char]
        else:
            # Raise a ValueError if a character is not found in the binary codes.
            raise ValueError(f"Character '{char}' not found in trie")

    return encoded_string  # Return the encoded string.

def print_node_list(water_network):
    for i in range(1, water_network._size + 1):
        print(water_network._node_list[i])
    

def print_edge_list(water_network):
    for i in range(1, water_network._size + 1):
        print(water_network._edge_list[i])
    



def main():
    water_network = read_csv_file("water_data.csv")
    visited = [False] * (water_network._size + 1)
    path = []
    # Q1
    print("\n\nQ1")
    cycles = []
    water_network.dfs_detect_cycles(2, visited, 2, path, cycles)
    print(cycles)
    # Q2
    print("\n\nQ2")
    prim = water_network.shortest_flight_path_search(Coordinates(200, 350), Coordinates(410, 600))
    print("Path:", " -> ".join(map(str, prim)))

    
    # Q3
    print("\n\nQ3")
    print(water_network.chemical_sources([(58,3),(55,10),(52,5)]))
    
    # Q4
    print("\n\nQ4")
    name_text = read_river_names_file("River_creek_names.csv")
    root = get_encoding_trie(name_text)
    encoded_string = code_string(root, "Daly")
    print(encoded_string)
    
    

if __name__ == "__main__":
    main()