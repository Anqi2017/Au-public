import string, sys, random

class Graph:
  # Use directed graph by default
  def __init__(self,directionless=False):
    self.__edges = {}
    self.__nodes = {}
    self.__directionless=directionless

  def get_edges(self):
    return    

  def get_nodes(self):
    return

  def add_node(self,node):
    self.__nodes[node.get_id()] = node
    return

  def get_children(self,node):
    if self.find_cycle() or self.__directionless:
      sys.stderr.write("ERROR: do cannot find a branch when there are cycles in the graph\n")
      sys.exit()
    v = self.__get_children(node.get_id())
    return [self.__nodes[i] for i in v]

  #assumes you have no cycles and its a directional graph
  def __get_children(self,nodeid):
    if nodeid not in self.__edges: #this is a leaf so traverse no farther
      return []
    kids = []
    for j in self.__edges[nodeid]:
      v = self.__get_children(j)
      for k in v:  
         if k not in kids: kids.append(k)
      kids.insert(0,j)
    return kids

  def get_roots(self):
    if self.__directionless:
      sys.stderr.write("ERROR: can't get roots of an undirected graph\n")
      sys.exit()
    outputids = self.__nodes.keys()
    #print outputids
    for j in self.__nodes:
      for ei in self.__edges:
        if j in self.__edges[ei]:
          outputids.remove(j)
          break
    return [self.__nodes[x] for x in outputids]

  def add_edge(self,edge,verbose=True):
    #make sure nodes are in the nodes
    if edge.get_node1().get_id() not in self.__nodes:
      self.__nodes[edge.get_node1().get_id()] = edge.get_node1()
    if edge.get_node2().get_id() not in self.__nodes:
      self.__nodes[edge.get_node2().get_id()] = edge.get_node2()
    # now add edge
    ids = edge.get_node_ids()    
    if ids[0] not in self.__edges:
      self.__edges[ids[0]] = set()
    if ids[1] in self.__edges[ids[0]] and verbose==True:
      sys.stderr.write("WARNING overwriting repeat edge.\n")
    self.__edges[ids[0]].add(ids[1])
    if edge.is_directionless():
      if ids[1] not in self.__edges:
        self.__ends[ids[1]] = set()
      if ids[0] in self.__edges[ids[1]] and verbose == True:
        sys.stderr.write("WARNING overwriting repeat edge.\n")
      self.__ends[ids[1]].add(ids[0])
    return

  def get_status_string(self):
    ostr = ''
    ostr += "----------------\n"
    ostr += "Node count: "+str(len(self.__nodes.keys()))+"\n"
    ostr += "Edge count: "+str(len(self.__edges.keys()))+"\n"
    if not self.__directionless:
      ostr += "Root count: "+str(len(self.get_roots()))+"\n"
    return ostr

  def remove_node(self,node):
    nid = node.get_id()
    #remove edges associated with this node
    toremove = []
    for i in self.__edges:
      for j in self.__edges[i]:
        if i == nid or j == nid:
          toremove.append([self.__nodes[i],self.__nodes[j]])
    for r in toremove:
      self.remove_edge(r[0],r[1])
    del self.__nodes[nid]

  # remove edge
  def remove_edge(self,node1,node2):
    if node1.get_id() in self.__edges:
      if node2.get_id() in self.__edges[node1.get_id()]:
        self.__edges[node1.get_id()].remove(node2.get_id())
      if len(self.__edges[node1.get_id()]) == 0:
        del self.__edges[node1.get_id()]
    if self.__directionless:
      if node2.get_id() in self.__edges: 
        if node1.get_id() in self.__edges[node2.get_id()]:
          self.__edges[node2.get_id()].remove(node1.get_id())
        if len(self.__edges[node2.get_id()]) == 0:
          del self.__edges[node2.get_id()]
    
  #remove cycles by mergine cyclic nodes into single nodes
  #their payloads are added to a list
  def merge_cycles(self):
    #delete any self cycles first
    for i in self.__edges:
      for j in self.__edges[i]:
        if self.__edges[i]==j: 
          self.remove_edge(self.__nodes[i],self.__nodes[j])
    while True:
      res = self.find_cycle()
      if not res:  return # we have finished.. there are no cycles
      # If we are here we need to merge
      if len(res) == 1: 
        sys.stderr.write("ERROR: Unexpected Self-cycle.\n")
        sys.exit()
      resids = [x.get_id() for x in res]
      # merge edges unless that edge is to one of the nodes we are removing
      for i in range(1,len(res)):
        for v in res[i].get_payload(): res[0].get_payload().append(v)
        if res[i].get_id() in self.__edges:
          for e2id in self.__edges[res[i].get_id()]:
            if e2id not in resids:
              self.add_edge(Edge(res[0],res[i]),verbose=False)
              if self.__directionless:
                self.add_edge(Edge(res[i],res[0]),verbose=False)
      # remove any nodes and edges connected to nodes we are removing
      for r in res[1:]:
        self.remove_node(r)

  #return a single cycle, greedy first one found
  #in terms of nodes return as an array of nodes or None
  def find_cycle(self):
    #depth first search through nodes
    for nid in self.__nodes.keys():
      res = self.__find_cycle_node([],nid)
      if res: return [self.__nodes[x] for x in res]
    return None

  #Internal function
  # Return the first cycle found form a starting node in terms of an
  # array of node ids          
  def __find_cycle_node(self,starting_ids,current):
    if current in starting_ids: return starting_ids
    if current not in self.__edges: return None
    for i in self.__edges[current]:
      newstarts = starting_ids[:]
      newstarts.append(current)
      res = self.__find_cycle_node(newstarts,i)
      if res: return res
    return None

# directed graph by default
class Edge:
  def __init__(self,node1,node2,directionless=False,weight=None):
    self.__node1 = node1
    self.__node2 = node2
    self.__directionless = directionless
    self.__weight = weight
  def set_weight(self,weight):  self.__weight = weight
  def get_weight(self): return self.__weight
  def get_node_ids(self): return [self.__node1.get_id(),self.__node2.get_id()]
  def is_directionless(self): return self.__directionless
  def get_node1(self): return self.__node1
  def get_node2(self): return self.__node2

# payload is a list of entries
class Node:
  def __init__(self,payload=None):
    self.__id = ''.join([random.choice(string.lowercase+string.uppercase+string.digits) for i in range(0,20)])
    self.__payload = []
    if payload:
      self.__payload = payload
  def get_payload(self):
    return self.__payload
  def set_payload(self,payload):
    self.__payload = payload
  def get_id(self):
    return self.__id
