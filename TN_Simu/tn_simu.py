import numpy as np 
from scipy import linalg
from qh_gates import get_matrix
from qh_circuit import QHCircuit, GateInstruction


class TNArchitecture:
    def __init__(self, nnodes):
        self.nnodes = nnodes
        self.edges = []
        self.neighbours = []
    
    @property
    def nedges(self):
        return len(self.edges)
    
    def compute_neighbours(self):
        self.neighbours = [[] for _ in range(self.nnodes)]
        for (i,j) in self.edges:
            self.neighbours[i].append(j)
            self.neighbours[j].append(i)
        self.neighbours = [tuple(each) for each in self.neighbours]


class PEPSArchitecture(TNArchitecture):
    def __init__(self, length, width):
        super().__init__(length * width)
        self.length = length # horizontal 
        self.width = width # vertical
        for x in range(length-1):
            for y in range(width-1):
                self.edges.append((length*y + x, length*y + x + 1))
                self.edges.append((length*y + x, length*y + x + length))
        for x in range(length-1):
            self.edges.append((length*(width-1) + x, length*(width-1) + x + 1))
        for y in range(width-1):
            self.edges.append((length*(y+1) - 1, length*(y+1) - 1 + length))
        self.compute_neighbours()


class MPSArchitecture(TNArchitecture):
    def __init__(self, length):
        super().__init__(length)
        self.length = length 
        for x in range(length-1):
            self.edges.append((x, x+1))
        self.compute_neighbours()


class IBM65Architecture(TNArchitecture):
    def __init__(self):
        super().__init__(65)
        self.edges = self.ibm65_edges()
        self.compute_neighbours

    def ibm65_edges(self):
        edges = [] 
        edges += [(i,i+1) for i in range(0, 9)]
        edges += [(i,i+1) for i in range(13, 23)]
        edges += [(i,i+1) for i in range(27, 37)]
        edges += [(i,i+1) for i in range(41, 51)]
        edges += [(i,i+1) for i in range(55, 64)]
        edges += [(0,10), (10,13), (4,11), (11,17), (8,12), (12,21)]
        edges += [(15,24), (24,29), (19,25), (25,33), (23,26), (26,37)]
        edges += [(27,38), (38,41), (31,39), (39,45), (35,40), (40,49)]
        edges += [(43,52), (52,56), (47,53), (53,60), (51,54), (54,64)]
        return edges


class IBM27Architecture(TNArchitecture):
    def __init__(self):
        super().__init__(27)
        self.edges = [(0,1),(1,2),(2,3),(1,4),(3,5),(4,7),(5,8),(6,7),(8,9),(7,10),(8,11),(10,12),(11,14),(12,13),(13,14),
                      (12,15),(14,16),(15,18),(16,19),(17,18),(19,20),(18,21),(19,22),(21,23),(22,25),(23,24),(24,25),(25,26)]
        self.compute_neighbours()


class IBM16Architecture(TNArchitecture):
    def __init__(self):
        super().__init__(16)
        self.edges = [(0,1),(1,2),(2,3),(1,4),(3,5),(4,7),(5,8),(6,7),(8,9),(7,10),(8,11),(10,12),(11,14),(12,13),(13,14),(12,15)]
        self.compute_neighbours()


class TNSimulator():
    """ 
    TN simulator takes in an architecture and a circuit. The circuit is assumed to have been 
    transpiled wrt to the architecture.

    Parameters:
        gamma: tensors.
        TODO: the followings are for tmp uses.
        classical_reg: classical register
        contraction_order: either 'default' or a list of integers representing the contraction order. 
            the contraction order for 'default' is list(range(nqubits))
    """
    def __init__(self, architecture: TNArchitecture, circuit: QHCircuit, xi=1e-4, max_chi=2**10):
        self.architecture = architecture
        self.neighbours = architecture.neighbours
        self.nqubits = architecture.nnodes
        self.circuit = circuit
        self.xi = xi 
        self.max_chi = max_chi
        self.verbose = 0
        self.gamma: list[np.ndarray] = []
        self.classical_reg = []
        self.contraction_order = 'default'
        self.initialize()


    def initialize(self):
        self.gamma = []
        for neib in self.neighbours:
            self.gamma.append(np.array([1,0], dtype=np.complex128).reshape((2,)+(1,)*len(neib)))

    
    def get_statistics(self):
        """
        Get the expected contraction FLOP.
        TODO: consider two-layer contraction.
        """
        def get_size(npara):
            b = 16 * npara # 128 bits per parameter
            if b < 1024:
                return f"{b} B"
            elif b < 1024 ** 2:
                return f"{b/1024:.0f} KiB"
            elif b < 1024 ** 3:
                return f"{b/1024**2:.0f} MiB"
            elif b < 1024 ** 4:
                return f"{b/1024**3:.0f} GiB"
            else:
                return f"{b/1024**4:.0f} TiB"
            
        degree_list = np.array([len(each) for each in self.neighbours])
        max_degree = degree_list.max()
        max_chi = self.max_chi
        max_n_para = 2 * (max_chi ** degree_list).sum()
        current_n_para = sum([np.prod(each.shape) for each in self.gamma])

        print(f"Highest degree: {max_degree}; Max bond dimension: {max_chi}")
        print(f"Max number of parameters:     {max_n_para} ({get_size(max_n_para)})")
        print(f"Current number of parameters: {current_n_para} ({get_size(current_n_para)})")


    def compute_contraction_flop(self):
        def locate_target(index_edge, target_node): 
            results = []
            for i in range(len(index_edge)):
                if index_edge[i][1] == target_node:
                    results.append((index_edge[i][0], i))
            if len(results) == 0:
                raise ValueError(f"Cannot find target node {target_node}")
            return results # a list of (source_node, idx)
        
        def prod(*arr):
            rst = 1
            for each in arr:
                rst *= each 
            return rst

        order = list(range(self.nqubits)) if self.contraction_order == 'default' else self.contraction_order

        first_node = order[0]
        shape_list = [each.shape[1:] for each in self.gamma]
        flops = 0
        tmp_dim = prod(*(shape_list[first_node]))
        tmp_edge_index = [(first_node, each) for each in self.neighbours[first_node]]
        for qubit in order[1:]:
            source_and_leg_idx = locate_target(tmp_edge_index, qubit)
            qubit_neighb = self.neighbours[qubit]
            qubit_shape = shape_list[qubit]
            qubit_dim = prod(*qubit_shape)
            shared_shape = [qubit_shape[qubit_neighb.index(source_node)] for source_node, _ in source_and_leg_idx]
            shared_dim = prod(*shared_shape)
            flops += tmp_dim * (qubit_dim // shared_dim) * (6 * shared_dim - 2)
            tmp_dim *= int(qubit_dim / shared_dim**2)

            tmp_edge_index += [(qubit, neib) for neib in self.neighbours[qubit]]
            for source, _ in source_and_leg_idx:
                tmp_edge_index.remove((source, qubit))
                tmp_edge_index.remove((qubit, source))

        return flops


    def get_amplitude(self, id: str|int):
        def locate_target(index_edge, target_node): 
            results = []
            for i in range(len(index_edge)):
                if index_edge[i][1] == target_node:
                    results.append((index_edge[i][0], i))
            if len(results) == 0:
                raise ValueError(f"Cannot find target node {index_edge} {target_node}")
            return results

        id = (r"{:0" + f"{self.nqubits}" + r"b}").format(id) if isinstance(id, int) else id
        id = [int(b) for b in id[::-1]]
        # print(id)
        
        order = list(range(self.nqubits)) if self.contraction_order == 'default' else self.contraction_order
        first_node = order[0]
        tmp = self.gamma[first_node][id[first_node]].copy()
        tmp_index_edge = [(first_node, each) for each in self.neighbours[0]]
        for next_node in order[1:]:
            next_neib = self.neighbours[next_node]
            n1 = len(tmp_index_edge)
            n2 = len(next_neib)
            tmp_idx = list(range(n1))
            next_idx = list(range(n1, n1+n2))

            tmp_index_edge += [(next_node, each) for each in next_neib]

            source_and_leg_idx = locate_target(tmp_index_edge, next_node)
            for count, (source_node, leg_idx) in enumerate(source_and_leg_idx):
                tmp_idx[leg_idx] = 51 - count # indices to be contracted starts from 51 reverse
                next_idx[next_neib.index(source_node)] = 51 - count
                tmp_index_edge.remove((source_node, next_node))
                tmp_index_edge.remove((next_node, source_node))

            n_contraction = len(source_and_leg_idx)
            out_idx = [each for each in tmp_idx if each < 51 - n_contraction] + [each for each in next_idx if each < 51 - n_contraction]

            tmp = np.einsum(tmp, tmp_idx, self.gamma[next_node][id[next_node]], next_idx, out_idx)

        return tmp


    def get_all_amplitudes(self, verbose=False):
        from tqdm import tqdm
        return np.array([self.get_amplitude(i) for i in tqdm(range(2**self.nqubits), disable=not verbose)])
    

    def get_probability(self, id: str|int):
        amp = self.get_amplitude(id)
        return np.abs(amp) ** 2


    def get_all_probabilities(self):
        amps = self.get_all_amplitudes()
        return np.abs(amps) ** 2


    def sample_amplitudes(self, nsamples=1e6, verbose=True):
        import random # np.random.randint cannot handle extremely large integers
        from tqdm import tqdm
        if 2**self.nqubits < nsamples:
            print(f'Warning: only {2**self.nqubits} are available.')
            return self.get_all_amplitudes()
        amps = []
        for idx in tqdm(range(int(nsamples)), disable=not verbose):
            idx = random.randint(0, 2**self.nqubits)
            amps.append(self.get_amplitude(idx))
        return np.array(amps)


    def single_qubit_gate(self, k, U):
        self.gamma[k] = np.einsum('im, m... -> i...', U, self.gamma[k])
    

    def get_measurement_probability(self, qubit, state=0):
        """
        TODO: re-write this code.
        Remark: This method is an early version. Rarely used later.
        """
        def locate_target(index_edge, target_node): 
            results = []
            for i in range(len(index_edge)):
                if index_edge[i][1] == target_node:
                    results.append((index_edge[i][0], i))
            if len(results) == 0:
                raise ValueError("Cannot find target node")
            return results

        gamma_norm = []
        for i, each in enumerate(self.gamma):
            shape_each = each.shape[1:] # ignore the physical dimension
            virtual_rank = len(shape_each)
            gamma_norm_shape = [int(shape_i ** 2) for shape_i in shape_each]
            no_conj_idx = list(range(virtual_rank))
            conj_idx = list(range(51, 51-virtual_rank, -1))    
            out_idx = sum(zip(no_conj_idx, conj_idx), start=())
            if i == qubit:
                tmp = np.einsum(each[state], no_conj_idx, each[state].conj(), conj_idx, out_idx, optimize=True)
                gamma_norm.append(tmp.reshape(gamma_norm_shape))
            else:
                tmp = np.einsum(each, [virtual_rank] + no_conj_idx, each.conj(), [virtual_rank] + conj_idx, out_idx, optimize=True)
                gamma_norm.append(tmp.reshape(gamma_norm_shape))

        order = list(range(self.nqubits)) if self.contraction_order == 'default' else self.contraction_order
        first_node = order[0]

        tmp = gamma_norm[first_node].copy()
        tmp_index_edge = [(first_node, each) for each in self.neighbours[first_node]]
        for next_node in order[1:]:
            next_neib = self.neighbours[next_node]
            n1 = len(tmp_index_edge)
            n2 = len(next_neib)
            tmp_idx = list(range(n1))
            next_idx = list(range(n1, n1+n2))

            tmp_index_edge += [(next_node, each) for each in next_neib]

            source_and_leg_idx = locate_target(tmp_index_edge, next_node)
            for count, (source_node, leg_idx) in enumerate(source_and_leg_idx):
                tmp_idx[leg_idx] = 51 - count # indices to be contracted starts from 51 reverse
                next_idx[next_neib.index(source_node)] = 51 - count
                tmp_index_edge.remove((source_node, next_node))
                tmp_index_edge.remove((next_node, source_node))

            n_contraction = len(source_and_leg_idx)
            out_idx = [each for each in tmp_idx if each < 51 - n_contraction] + [each for each in next_idx if each < 51 - n_contraction]

            tmp = np.einsum(tmp, tmp_idx, gamma_norm[next_node], next_idx, out_idx, optimize=True)

        return np.abs(tmp)


    def sample_bitstring(self, repetition, verbose=False):
        """
        Return: list<string>
            Each string is a length-n binary string representing measurement outcomes.
        """
        from tqdm import tqdm
        bitstrings = []
        for _ in tqdm(range(repetition), disable=not verbose):
            original = [each.copy() for each in self.gamma]
            bitstr = ''
            for qubit in range(self.nqubits):
                prob = self.get_measurement_probability(qubit, state=0)
                output = int(np.random.rand() > prob)
                if output == 0:
                    self.single_qubit_gate(qubit, np.array([[prob**-0.5, 0], [0, 0]]))
                else:
                    self.single_qubit_gate(qubit, np.array([[0, 0], [0, (1-prob)**-0.5]]))
                bitstr = str(output) + bitstr
            bitstrings.append(bitstr)
            self.gamma = original

        return bitstrings
        

    def two_qubit_gate_merge_split(self, k, l, U):
        try:
            idx_l_in_k = self.neighbours[k].index(l)
            idx_k_in_l = self.neighbours[l].index(k)
        except:
            raise ValueError(f"Qubits {k} and {l} are not connected!")
        
        def stage_1_transpose_key(dimension, bond_idx):
            return list(range(1, 1+bond_idx)) + list(range(2+bond_idx, dimension)) + [0, 1+bond_idx]   

        def stage_2_transpose_key(dimension, bond_idx):
            return list(range(bond_idx+1)) + [dimension-1] + list(range(bond_idx+1, dimension-1))

        shape_k = list(self.gamma[k].shape)
        shape_l = list(self.gamma[l].shape)

        dim_k = len(shape_k)
        dim_l = len(shape_l)

        total_dim_k = np.prod(shape_k)
        total_dim_l = np.prod(shape_l)

        chi_local = shape_k[idx_l_in_k+1]
        
        # Eeshape 
        gamma_k = np.transpose(self.gamma[k], stage_1_transpose_key(dim_k, idx_l_in_k)).reshape(-1, 2, chi_local) # (k1 k2 k3), (m), (alpha)
        gamma_l = np.transpose(self.gamma[l], stage_1_transpose_key(dim_l, idx_k_in_l)).reshape(-1, 2, chi_local) # (l1 l2 l3), (n), (alpha)

        Theta = np.einsum('ijmn, kma, lna -> ikjl', U, gamma_k, gamma_l).reshape(total_dim_k//chi_local, -1) # (i k1 k2 k3) (j l1 l2 l3)

        # SVD
        u, s, vh = np.linalg.svd(Theta)
        propose_new_chi_local = (s > self.xi).sum()
        new_chi_local = min(propose_new_chi_local, self.max_chi)

        s *= ((s ** 2).sum() / (s[:new_chi_local] ** 2).sum() ) ** 0.5 # rescale

        # Truncate, calculate new dimension
        updated_gamma_k = (u[:, :new_chi_local] * s[None, :new_chi_local]).reshape(2, -1, new_chi_local) # (i) (k1 k2 k3) (t)

        updated_gamma_l = vh[:new_chi_local, :].reshape(new_chi_local, 2, -1) # (t) (j) (l1 l2 l3)
        updated_gamma_l = updated_gamma_l.transpose((1,2,0)) # (j) (l1 l2 l3) (t)

        shape_k_b4_state2_tranpose = shape_k[:idx_l_in_k+1] + shape_k[idx_l_in_k+2:] + [new_chi_local]
        shape_l_b4_state2_tranpose = shape_l[:idx_k_in_l+1] + shape_l[idx_k_in_l+2:] + [new_chi_local]

        updated_gamma_k = updated_gamma_k.reshape(shape_k_b4_state2_tranpose)
        updated_gamma_l = updated_gamma_l.reshape(shape_l_b4_state2_tranpose)

        self.gamma[k] = np.transpose(updated_gamma_k, stage_2_transpose_key(dim_k, idx_l_in_k))
        self.gamma[l] = np.transpose(updated_gamma_l, stage_2_transpose_key(dim_l, idx_k_in_l))


    def two_qubit_gate_qrsvd(self, k, l, U):
        try:
            idx_l_in_k = self.neighbours[k].index(l)
            idx_k_in_l = self.neighbours[l].index(k)
        except:
            raise ValueError(f"Qubits {k} and {l} are not connected!")
        
        def stage_1_transpose_key(dimension, bond_idx):
            return list(range(1, 1+bond_idx)) + list(range(2+bond_idx, dimension)) + [0, 1+bond_idx]
        
        def stage_2_transpose_key(dimension, bond_idx):
            return list(range(bond_idx+1)) + [dimension-1] + list(range(bond_idx+1, dimension-1))

        shape_k = list(self.gamma[k].shape)
        shape_l = list(self.gamma[l].shape)

        dim_k = len(shape_k)
        dim_l = len(shape_l)

        total_dim_k = np.prod(shape_k)
        total_dim_l = np.prod(shape_l)

        chi_local = shape_k[idx_l_in_k+1]

        total_dim_k_no_l = total_dim_k // chi_local // 2
        total_dim_l_no_k = total_dim_l // chi_local // 2
        
        gamma_k = np.transpose(self.gamma[k], stage_1_transpose_key(dim_k, idx_l_in_k)).reshape(-1, 2*chi_local) # (k1 k2 k3), (m alpha)
        gamma_l = np.transpose(self.gamma[l], stage_1_transpose_key(dim_l, idx_k_in_l)).reshape(-1, 2*chi_local) # (l1 l2 l3), (n alpha)

        q_k, r_k = linalg.rq(gamma_k) # q_k: (k1 k2 k3), (beta);  r_k: (beta),  (m alpha)
        q_l, r_l = linalg.rq(gamma_l) # q_l: (l1 l2 l3), (gamma); r_l: (gamma), (n alpha)

        r_k = r_k.reshape(2*chi_local, 2, chi_local)
        r_l = r_l.reshape(2*chi_local, 2, chi_local)
        Theta = np.einsum('ijmn, bma, cna -> ibjc', U, r_k, r_l).reshape(4*chi_local, 4*chi_local) # (i beta) (j gamma)

        u, s, vh = np.linalg.svd(Theta) # u: (i beta) (t); vh: (t) (j gamma)
        propose_new_chi_local = (s > self.xi).sum()
        new_chi_local = min(propose_new_chi_local, self.max_chi)

        # Theta1 = Theta.reshape(2, 2*chi_local, 2, 2*chi_local)

        s *= ((s ** 2).sum() / (s[:new_chi_local] ** 2).sum() ) ** 0.5

        updated_r_k = (u[:, :new_chi_local] * s[None, :new_chi_local]).reshape(2, 2*chi_local, new_chi_local) # (i) (beta) (t)
        updated_r_l = vh[:new_chi_local, :].reshape(new_chi_local, 2, 2*chi_local) # (t) (j) (gamma)

        # Theta2 = np.einsum('ibt, tjc -> ibjc', updated_r_k, updated_r_l)

        shape_k_b4_state2_tranpose = shape_k[:idx_l_in_k+1] + shape_k[idx_l_in_k+2:] + [new_chi_local]
        shape_l_b4_state2_tranpose = shape_l[:idx_k_in_l+1] + shape_l[idx_k_in_l+2:] + [new_chi_local]
        updated_gamma_k = np.einsum('kb, ibt -> ikt', q_k.reshape(total_dim_k_no_l, 2*chi_local), updated_r_k).reshape(shape_k_b4_state2_tranpose)
        updated_gamma_l = np.einsum('kb, tib -> ikt', q_l.reshape(total_dim_l_no_k, 2*chi_local), updated_r_l).reshape(shape_l_b4_state2_tranpose)
        
        self.gamma[k] = np.transpose(updated_gamma_k, stage_2_transpose_key(dim_k, idx_l_in_k))
        self.gamma[l] = np.transpose(updated_gamma_l, stage_2_transpose_key(dim_l, idx_k_in_l))


    def two_qubit_gate_decomposed(self, k, l, U_decomp):
        raise NotImplementedError()


    def apply_instruction(self, ins: GateInstruction, method='qr-svd'):
        if ins.type == GateInstruction.UNITARY_GATE:
            if ins.nqubits == 1:
                self.single_qubit_gate(ins.t_qubit[0], get_matrix(ins.gate, *(ins.params)))
            elif ins.nqubits == 2:
                more_sig, less_sig = ins.t_qubit + ins.c_qubit
                mat = get_matrix(ins.gate, *(ins.params))
                if method == 'qr-svd':
                    self.two_qubit_gate_qrsvd(more_sig, less_sig, mat)
                elif method == 'merge-split':
                    self.two_qubit_gate_merge_split(more_sig, less_sig, mat)
                else:
                    raise NotImplementedError(f"Unknown method {method}. Supported are 'qr-svd' and 'merge-split'.")
            else:
                print(f"Warning: cannot apply instrution {ins.__str__()}")
        elif ins.type == GateInstruction.MEASUREMENT:
            qubit = ins.t_qubit[0]
            prob0 = self.get_measurement_probability(qubit, state=0)
            if np.random.rand() < prob0:
                # measured 0
                self.classical_reg.append(0)
                self.single_qubit_gate(qubit, np.array([[prob0**-0.5, 0], [0, 0]]))
            else:
                # measured 1 
                self.classical_reg.append(1)
                self.single_qubit_gate(qubit, np.array([[0, 0], [0, (1-prob0)**-0.5]]))
        elif ins.type == GateInstruction.COLLAPSE:
            qubit = ins.t_qubit[0]
            state = ins.params[0]
            self.collapse(qubit, state)


    def collapse(self, qubit, state):
        prob0 = self.get_measurement_probability(qubit, state=0)
        prob1 = 1 - prob0
        if state == 0:
            mat = np.array([[prob0**-0.5, 0], [0, 0]]) if prob0 > prob1 else np.array([[0, prob1**-0.5], [0, 0]])
        else: # state == 1
            mat = np.array([[0, 0], [0, prob1**-0.5]]) if prob1 > prob0 else np.array([[0, 0], [prob0**-0.5, 0]])
        # print(mat)
        self.single_qubit_gate(qubit, mat)


    def simulate(self, method='qr-svd', verbose=False):
        """
        Input:
            method: currently support 'qr-svd' and 'merge-split'
        """
        from tqdm import tqdm
        for ins in tqdm(self.circuit.instructions, disable=not verbose):
            self.apply_instruction(ins, method=method)