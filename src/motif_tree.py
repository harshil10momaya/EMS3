import sys
from typing import Dict, Any, List, Optional
from abc import ABC, abstractmethod
from utils import WILDCARD_CODE, DNA_DOMAIN

# Define node allocators based on the C++ structures
def allocate_node_fast() -> Dict[str, Any]:
    """Node for MotifTreeFast: array-based children, heavy sharing info."""
    return {'children': [None] * len(DNA_DOMAIN), 'sharing_info': 0}

def allocate_node_slow() -> Dict[str, Any]:
    """Node for MotifTreeSlow equivalent: list-based children."""
    return {'children': [], 'sharing_info': 0}


# ====================================================================
# BASE MOTIF TREE (CRTP equivalent)
# ====================================================================

class MotifTreeBase(ABC):
    def __init__(self, max_depth: int, motifs: List[str], name: str, alloc_func: callable):
        self.domain: str = DNA_DOMAIN
        self.domain_size: int = len(DNA_DOMAIN)
        self.mask: int = (1 << self.domain_size) - 1 # 0b1111
        self.node_allocator = alloc_func
        self.root: Dict[str, Any] = self.node_allocator()
        self.max_depth: int = max_depth
        self.x: List[str] = [''] * max_depth # Temporary array for building motifs
        self.name: str = name
        self.motifs: List[str] = motifs

    @abstractmethod
    def traverse_recursive(self, node: Dict[str, Any], depth: int): pass
    def traverse(self):
        self.motifs.clear()
        self.traverse_recursive(self.root, 0)

    @abstractmethod
    def insert_recursive(self, node: Dict[str, Any], motif: str, depth: int): pass
    def insert(self, motif: str): self.insert_recursive(self.root, motif, 0)

    @abstractmethod
    def intersect_recursive(self, to_node: Dict[str, Any], from_node: Dict[str, Any], depth: int) -> Optional[Dict[str, Any]]: pass
    def intersect(self, other: 'MotifTreeBase'):
        new_root = self.intersect_recursive(self.root, other.root, 0)
        self.root = new_root if new_root else self.node_allocator()


# ====================================================================
# MOTIF TREE FAST (C++ `MotifTreeFast` equivalent)
# ====================================================================

class MotifTreeFast(MotifTreeBase):
    def __init__(self, max_depth: int, motifs: List[str], name: str):
        super().__init__(max_depth, motifs, name, allocate_node_fast)

    def _empty_node(self, node: Dict[str, Any]) -> bool:
        return all(child is None for child in node['children'])
    
    def traverse_recursive(self, node: Dict[str, Any], depth: int):
        if depth == self.max_depth:
            self.motifs.append("".join(self.x))
            return

        for i in range(self.domain_size):
            child_node = node['children'][i]
            if child_node is not None:
                self.x[depth] = str(i) 
                self.traverse_recursive(child_node, depth + 1)

    def insert_recursive(self, node: Dict[str, Any], motif: str, depth: int):
        if depth >= self.max_depth: return
        ch = int(motif[depth]) if motif[depth].isdigit() else WILDCARD_CODE

        if ch == WILDCARD_CODE:
            current_info = 0
            for j in range(self.domain_size):
                child_node = node['children'][j]
                if child_node is not None and not (current_info & (1 << j)):
                    self.insert_recursive(child_node, motif, depth + 1)
                    current_info |= child_node['sharing_info'] 
            
            remaining = (~current_info) & self.mask
            if remaining:
                new_node = self.node_allocator()
                new_node['sharing_info'] = remaining
                for j in range(self.domain_size):
                    if (remaining & (1 << j)): node['children'][j] = new_node
                self.insert_recursive(new_node, motif, depth + 1)
        else:
            child_mask = (1 << ch)
            child_node = node['children'][ch]
            
            if child_node is None:
                new_node = self.node_allocator()
                new_node['sharing_info'] = child_mask
                node['children'][ch] = new_node
                self.insert_recursive(new_node, motif, depth + 1)
            else:
                old_info = child_node['sharing_info']
                if old_info & (~child_mask):
                    old_node = child_node 
                    new_node = self.node_allocator()
                    new_node['sharing_info'] = child_mask
                    old_node['sharing_info'] &= ~child_mask 
                    node['children'][ch] = new_node 
                    self.insert_recursive(new_node, motif, depth + 1)
                else: 
                    self.insert_recursive(child_node, motif, depth + 1)

    def intersect_recursive(self, to_node: Dict[str, Any], from_node: Optional[Dict[str, Any]], depth: int) -> Optional[Dict[str, Any]]:
        if depth >= self.max_depth: return to_node 
        if from_node is None: return None

        new_to_node = self.node_allocator()

        for j in range(self.domain_size):
            to_child_original = to_node['children'][j]
            from_child_original = from_node['children'][j]
            
            if to_child_original is None: continue

            to_info_original = to_child_original['sharing_info']
            from_info = from_child_original['sharing_info'] if from_child_original else 0
            common = to_info_original & from_info

            if common:
                intersected_child = self.intersect_recursive(to_child_original, from_child_original, depth + 1)

                if intersected_child is not None:
                    intersected_child['sharing_info'] = common
                    for k in range(self.domain_size):
                        if common & (1 << k):
                            new_to_node['children'][k] = intersected_child
                            new_to_node['sharing_info'] |= (1 << k)

        if self._empty_node(new_to_node) and depth < self.max_depth: return None
        return new_to_node

# ====================================================================
# MOTIF TREE SIMPLE (C++ `MotifTreeSlow` equivalent for ems2m)
# ====================================================================

class MotifTreeSimple(MotifTreeBase):
    def __init__(self, max_depth: int, motifs: List[str], name: str):
        super().__init__(max_depth, motifs, name, allocate_node_slow)
    
    def traverse_recursive(self, node: Dict[str, Any], depth: int):
        if depth == self.max_depth:
            self.motifs.append("".join(self.x))
            return

        children_map = {}
        for child in node['children']:
            info = child['sharing_info']
            for j in range(self.domain_size):
                if info & (1 << j): children_map[j] = child
                    
        for i in range(self.domain_size):
            child_node = children_map.get(i)
            if child_node is not None:
                self.x[depth] = str(i) 
                self.traverse_recursive(child_node, depth + 1)

    def insert_recursive(self, node: Dict[str, Any], motif: str, depth: int):
        if depth >= self.max_depth: return

        ch = int(motif[depth]) if motif[depth].isdigit() else WILDCARD_CODE
        target_info = self.mask if ch == WILDCARD_CODE else (1 << ch)
        
        found_child = None
        for child in node['children']:
            if child['sharing_info'] & target_info:
                found_child = child
                break
                
        if found_child:
            self.insert_recursive(found_child, motif, depth + 1)
        else:
            new_node = self.node_allocator()
            new_node['sharing_info'] = target_info
            node['children'].append(new_node)
            self.insert_recursive(new_node, motif, depth + 1)

    def intersect_recursive(self, to_node: Dict[str, Any], from_node: Dict[str, Any], depth: int) -> Optional[Dict[str, Any]]:
        if depth >= self.max_depth: return to_node 
        if not from_node['children']: return None

        new_children = []
        
        for to_child in to_node['children']:
            for from_child in from_node['children']:
                common = to_child['sharing_info'] & from_child['sharing_info']
                
                if common:
                    new_node = self.node_allocator()
                    intersected_child = self.intersect_recursive(to_child, from_child, depth + 1)
                    
                    if intersected_child:
                        new_node['sharing_info'] = common
                        new_node['children'] = intersected_child['children']
                        new_children.append(new_node)

        if not new_children and depth < self.max_depth: return None
        
        new_to_node = self.node_allocator()
        new_to_node['children'] = new_children
        return new_to_node