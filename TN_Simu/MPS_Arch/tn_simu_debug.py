import numpy as np 
from qh_gates import get_matrix
from qh_circuit import QHCircuit, GateInstruction


class MPS_Simulator:
    def __init__(self, circuit: QHCircuit, xi=1e-4, max_chi=None):
        self.verbose = False
        self.circuit = circuit
        self.n = circuit.nqubits
        self.xi = xi
        self.max_chi = 2 ** (self.n // 2) if max_chi is None else max_chi
        self.initialize()
        self.classical_reg = []

    def initialize(self):
        self._gamma = [np.zeros((2, 1, 1), dtype=np.complex128) for _ in range(self.n)]
        for gamma in self._gamma:
            gamma[0, 0, 0] = 1 


    def get_amplitude(self, id: str|int):
        index = int(id, 2) if isinstance(id, str) else id

        tmp = self._gamma[0][index%2].copy()
        # print(tmp.shape)
        index //= 2
        for i in range(1, self.n):
            tmp = tmp @ self._gamma[i][index%2]
            index //= 2

        return tmp[0,0]

    def get_all_amplitudes(self):
        return np.array([self.get_amplitude(i) for i in range(2**self.n)])
    
    def get_probability(self, id: str|int):
        amp = self.get_amplitude(id)
        return np.abs(amp) ** 2

    def get_all_probabilities(self):
        amps = self.get_all_amplitudes()
        return np.abs(amps) ** 2

    def get_measurement_probability(self, qubit, state=0):
        tmp = np.array([1], dtype=np.complex128).reshape((1,1))
        for i in range(self.n):
            if i == qubit:
                tmp = np.einsum('ab, ac, bd -> cd', tmp, self._gamma[i][state], self._gamma[i][state].conj())
            else:
                tmp = np.einsum('ab, iac, ibd -> cd', tmp, self._gamma[i], self._gamma[i].conj())
        
        return np.abs(tmp[0,0])


    def single_qubit_unitary(self, k, U):
        self._gamma[k] = np.einsum('im, mab -> iab', U, self._gamma[k])

    def two_qubit_unitary_consecutive(self, k, U):
        # k is more significant
        l = k + 1
        chi_left = self._gamma[k].shape[1]
        chi_local = self._gamma[k].shape[2]
        chi_right = self._gamma[l].shape[2]
        
        tmp = np.einsum('ijmn, mab, nbc -> iajc', U, self._gamma[k], self._gamma[l]) 
        tmp = tmp.reshape((2*chi_left, 2*chi_right)) 
        u, s, vh = np.linalg.svd(tmp)

        # check threshold
        thres_chi = (s > self.xi).sum()
        new_chi = min(thres_chi, self.max_chi)
        if new_chi > chi_local:
            chi_local = new_chi
            if self.verbose:
                print(f"Bond dimension increased to {chi_local}")

        s *= s.sum() / s[:chi_local].sum()

        self._gamma[k] = (u[:,:chi_local] * s[None, :chi_local]).reshape((2, chi_left, chi_local))
        self._gamma[l] = np.einsum('bic -> ibc', vh[:chi_local,:].reshape((chi_local, 2, chi_right)))


    def two_qubit_unitary(self, k, l, U: np.ndarray):
        if k+1 == l:
            self.two_qubit_unitary_consecutive(k, U)
        elif k == l+1:
            self.two_qubit_unitary_consecutive(l, np.einsum('ijmn -> jinm', U))
        else:
            raise ValueError(f"Qubit {k} and {l} are not connected!")
    

    def apply_instruction(self, ins: GateInstruction):
        if ins.type == GateInstruction.UNITARY_GATE:
            if ins.nqubits == 1:
                self.single_qubit_unitary(ins.t_qubit[0], get_matrix(ins.gate, *(ins.params)))
            elif ins.nqubits == 2:
                more_sig, less_sig = ins.t_qubit + ins.c_qubit
                mat = get_matrix(ins.gate, *(ins.params))
                self.two_qubit_unitary(more_sig, less_sig, mat)
            else:
                print(f"Warning: cannot apply instrution {ins.__str__()}")
        elif ins.type == GateInstruction.MEASUREMENT:
            qubit = ins.t_qubit[0]
            prob0 = self.get_measurement_probability(qubit, state=0)
            if np.random.rand() < prob0:
                # measured 0
                self.classical_reg.append(0)
                self.single_qubit_unitary(qubit, np.array([[prob0**-0.5, 0], [0, 0]]))
            else:
                # measured 1 
                self.classical_reg.append(1)
                self.single_qubit_unitary(qubit, np.array([[0, 0], [0, (1-prob0)**-0.5]]))
        elif ins.type == GateInstruction.COLLAPSE:
            qubit = ins.t_qubit[0]
            state = ins.params[0]
            prob0 = self.get_measurement_probability(qubit, state=0)
            prob1 = 1 - prob0
            if state == 0:
                mat = np.array([[prob0**-0.5, 0], [0, 0]]) if prob0 > prob1 else np.array([[0, prob1**-0.5], [0, 0]])
            else: # state == 1
                mat = np.array([[0, 0], [0, prob1**-0.5]]) if prob1 > prob0 else np.array([[0, 0], [prob0**-0.5, 0]])
            
            self.single_qubit_unitary(qubit, mat)

    def simulate(self):
        for ins in self.circuit.instructions:
            self.apply_instruction(ins)

        


class TN_Architecture:
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


class PEPS_Architecture(TN_Architecture):
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







class TN_Simulator():
    def __init__(self, architecture: TN_Architecture, circuit: QHCircuit, xi=1e-4, max_chi=2**10):
        self.architecture = architecture
        self.neighbours = architecture.neighbours
        self.nqubits = architecture.nnodes
        self.circuit = circuit
        self.xi = xi 
        self.max_chi = max_chi
        self.verbose = 0
        self.gamma = []
        self.classical_reg = []
        self._initialize_gamma()

    def _initialize_gamma(self):
        for neib in self.architecture.neighbours:
            self.gamma.append(np.array([1, 0j]).reshape((2,) + (1,) * len(neib)))
    
    
    def get_amplitude(self, id: str|int):
        def locate_target(index_edge, target_node): 
            results = []
            for i in range(len(index_edge)):
                if index_edge[i][1] == target_node:
                    results.append((index_edge[i][0], i))
            if len(results) == 0:
                raise ValueError("Cannot find target node")
            return results

        id = int(id, 2) if isinstance(id, str) else id
        
        tmp = self.gamma[0][id%2].copy()
        tmp_index_edge = [(0, each) for each in self.neighbours[0]]
        id //= 2
        for next_node in range(1, self.nqubits):
            next_neib = self.neighbours[next_node]
            n1 = len(tmp_index_edge)
            n2 = len(next_neib)
            tmp_idx = list(range(n1))
            next_idx = list(range(n1, n1+n2))

            if self.verbose > 2:
                print(f"[Debug 2] next node = {next_node}, physical = {id%2}")
                print(f"[Debug 2] index edge {tmp_index_edge} next neib {next_neib}")

            tmp_index_edge += [(next_node, each) for each in next_neib]

            source_and_leg_idx = locate_target(tmp_index_edge, next_node)
            for count, (source_node, leg_idx) in enumerate(source_and_leg_idx):
                tmp_idx[leg_idx] = 51 - count # indices to be contracted starts from 51 reverse
                next_idx[next_neib.index(source_node)] = 51 - count
                tmp_index_edge.remove((source_node, next_node))
                tmp_index_edge.remove((next_node, source_node))

            n_contraction = len(source_and_leg_idx)
            out_idx = [each for each in tmp_idx if each < 51 - n_contraction] + [each for each in next_idx if each < 51 - n_contraction]

            if self.verbose > 2:
                print(f"[Debug 2] tmp idx {tmp_idx}, next idx {next_idx}, out idx {out_idx}")

            tmp = np.einsum(tmp, tmp_idx, self.gamma[next_node][id%2], next_idx, out_idx)
            id //= 2

            if self.verbose > 2:
                print(f"[Debug 2] Contracted {next_node}. tmp shape {tmp.shape}")
                print(f"[Debug 2]\t index edge {tmp_index_edge}")

        return tmp
    
    def get_measurement_probability(self, qubit, state=0):
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
                tmp = np.einsum(each[state], no_conj_idx, each[state].conj(), conj_idx, out_idx)
                gamma_norm.append(tmp.reshape(gamma_norm_shape))
            else:
                tmp = np.einsum(each, [virtual_rank] + no_conj_idx, each.conj(), [virtual_rank] + conj_idx, out_idx)
                gamma_norm.append(tmp.reshape(gamma_norm_shape))


        tmp = gamma_norm[0].copy()
        tmp_index_edge = [(0, each) for each in self.neighbours[0]]
        for next_node in range(1, self.nqubits):
            next_neib = self.neighbours[next_node]
            n1 = len(tmp_index_edge)
            n2 = len(next_neib)
            tmp_idx = list(range(n1))
            next_idx = list(range(n1, n1+n2))

            if self.verbose > 2:
                print(f"[Debug 2] next node = {next_node}")
                print(f"[Debug 2] index edge {tmp_index_edge} next neib {next_neib}")

            tmp_index_edge += [(next_node, each) for each in next_neib]

            source_and_leg_idx = locate_target(tmp_index_edge, next_node)
            for count, (source_node, leg_idx) in enumerate(source_and_leg_idx):
                tmp_idx[leg_idx] = 51 - count # indices to be contracted starts from 51 reverse
                next_idx[next_neib.index(source_node)] = 51 - count
                tmp_index_edge.remove((source_node, next_node))
                tmp_index_edge.remove((next_node, source_node))

            n_contraction = len(source_and_leg_idx)
            out_idx = [each for each in tmp_idx if each < 51 - n_contraction] + [each for each in next_idx if each < 51 - n_contraction]

            if self.verbose > 2:
                print(f"[Debug 2] tmp idx {tmp_idx}, next idx {next_idx}, out idx {out_idx}")

            tmp = np.einsum(tmp, tmp_idx, gamma_norm[next_node], next_idx, out_idx)

            if self.verbose > 2:
                print(f"[Debug 2] Contracted {next_node}. tmp shape {tmp.shape}")
                print(f"[Debug 2]\t index edge {tmp_index_edge}")

        return np.abs(tmp)

    def get_all_amplitudes(self):
        return np.array([self.get_amplitude(i) for i in range(2**self.nqubits)])
    
    def get_probability(self, id: str|int):
        amp = self.get_amplitude(id)
        return np.abs(amp) ** 2

    def get_all_probabilities(self):
        amps = self.get_all_amplitudes()
        return np.abs(amps) ** 2


    def single_qubit_unitary(self, k, U):
        self.gamma[k] = np.einsum('im, m... -> i...', U, self.gamma[k])
    
    def two_qubit_unitary(self, k, l, U):
        def prod(*arr):
            p = 1
            for each in arr:
                p *= each 
            return p
        
        if not l in self.neighbours[k]:
            raise ValueError(f"Qubit {k} and {l} are not connected.")
    
        if self.verbose > 1:
            print(f"[Debug 1] Apply two-qubit gate on k={k} (more significant) and l={l} (less significant).")
    
        nk = len(self.neighbours[k])
        nl = len(self.neighbours[l])
        k_idx = list(range(nk)) 
        l_idx = list(range(nk, nk+nl)) 
        k_phy_idx_in = nk+nl
        l_phy_idx_in = nk+nl+1
        k_phy_idx_out = nk+nl+2
        l_phy_idx_out = nk+nl+3
        idx_l_in_k = self.neighbours[k].index(l)
        idx_k_in_l = self.neighbours[l].index(k)
        k_idx[idx_l_in_k] = 51
        l_idx[idx_k_in_l] = 51


        gate_idx = [k_phy_idx_out, l_phy_idx_out, k_phy_idx_in, l_phy_idx_in]
        out_idx = [k_phy_idx_out] + [each for each in k_idx if each != 51] + [l_phy_idx_out] + [each for each in l_idx if each != 51]

        k_idx = [k_phy_idx_in] + k_idx
        l_idx = [l_phy_idx_in] + l_idx

        shape_k = list(self.gamma[k].shape)
        shape_l = list(self.gamma[l].shape)

        if self.verbose > 2:
            tmp_dict = {}
            for i, neib in enumerate(self.neighbours[k]):
                tmp_dict[neib] = shape_k[i+1]
            print(f"[Debug 2] k={k} has {nk} neighbors. [qubit: chi] as follows: {tmp_dict}")

            tmp_dict = {}
            for i, neib in enumerate(self.neighbours[l]):
                tmp_dict[neib] = shape_l[i+1]
            print(f"[Debug 2] l={l} has {nl} neighbors. [qubit: chi] as follows: {tmp_dict}")

        shape_k_without_l = shape_k[:idx_l_in_k+1] + shape_k[idx_l_in_k+2:]
        shape_l_without_k = shape_l[:idx_k_in_l+1] + shape_l[idx_k_in_l+2:]
        dim_left = prod(*shape_k_without_l)
        dim_right = prod(*shape_l_without_k)
        dim_local = shape_k[idx_l_in_k+1]

        if self.verbose > 2:
            print(f"[Debug 2] Shape k without l: {shape_k_without_l}. Dimension Left: {dim_left}")
            print(f"[Debug 2] Shape l without k: {shape_l_without_k}. Dimension Right: {dim_right}")
            print(f"[Debug 2] Dimension Local {dim_local}")


        if self.verbose > 2:
            print(f"[Debug 2] Build Theta: k_idx {k_idx}, l_idx {l_idx}, gate idx {gate_idx}, out idx {out_idx}")
            print(f"[Debug 2] Among these, kl in {k_phy_idx_in, l_phy_idx_in}; kl out {k_phy_idx_out, l_phy_idx_out}")

        tmp = np.einsum(self.gamma[k], k_idx, self.gamma[l], l_idx, U, gate_idx, out_idx).reshape(dim_left, dim_right)

        # print(tmp)

        if self.verbose > 1:
            print(f"[Debug 1] Theta is reshaped to {tmp.shape}")


        u, s, vh = np.linalg.svd(tmp)
        
        if self.verbose > 2:
            print(f"[Debug 2] SVD output: U shape {u.shape}, V shape {vh.shape}")
            print(f"[Debug 2] Singular values: {s.round(3)}")

        thres_chi = (s > self.xi).sum()
        new_chi = min(thres_chi, self.max_chi)
        if new_chi > dim_local:
            dim_local = new_chi
            if self.verbose > 0:
                print(f"Bond dimension increased to {dim_local}")

        if self.verbose > 1:
            print(f"[Debug 1] Truncation level (local bond dimension): {dim_local}")

        s *= s.sum() / s[:dim_local].sum()

        if self.verbose > 2:
            print(f"[Debug 2] Truncated Singular value after normalise: {s[:dim_local]}")

        new_shape_k = shape_k[:idx_l_in_k+1] + shape_k[idx_l_in_k+2:] + [dim_local]

        new_gamma_k = (u[:,:dim_local] * s[None, :dim_local]).reshape(new_shape_k)

        idx_shuffle_k_in = list(range(nk)) + [51]
        idx_shuffle_k_out = list(range(nk))
        idx_shuffle_k_out.insert(idx_l_in_k+1, 51)

        self.gamma[k] = np.einsum(new_gamma_k, idx_shuffle_k_in, idx_shuffle_k_out)

        if self.verbose > 2:
            print(f"[Debug 2] gamma[k] transpose key {idx_shuffle_k_in} -> {idx_shuffle_k_out}")
            print(f"[Debug 2] gamma[k] dimension {shape_k} -> {new_shape_k} -> {self.gamma[k].shape}")

        new_shape_l = [dim_local] + shape_l[:idx_k_in_l+1] + shape_l[idx_k_in_l+2:]

        if self.verbose > 2:
            print(f"[Debug 2] gamma[l] dimension before tranpose {shape_l} -> {new_shape_l}")

        new_gamma_l = vh[:dim_local, :].reshape(new_shape_l)

        idx_shuffle_l_in = [51] + list(range(nl))
        idx_shuffle_l_out = list(range(nl))
        idx_shuffle_l_out.insert(idx_k_in_l+1, 51)

        if self.verbose > 2:
            print(f"[Debug 2] gamma[l] transpose key {idx_shuffle_l_in} -> {idx_shuffle_l_out}")

        self.gamma[l] = np.einsum(new_gamma_l, idx_shuffle_l_in, idx_shuffle_l_out)
    
    def apply_instruction(self, ins: GateInstruction):
        if ins.type == GateInstruction.UNITARY_GATE:
            if ins.nqubits == 1:
                self.single_qubit_unitary(ins.t_qubit[0], get_matrix(ins.gate, *(ins.params)))
            elif ins.nqubits == 2:
                more_sig, less_sig = ins.t_qubit + ins.c_qubit
                mat = get_matrix(ins.gate, *(ins.params))
                self.two_qubit_unitary(more_sig, less_sig, mat)
            else:
                print(f"Warning: cannot apply instrution {ins.__str__()}")
        elif ins.type == GateInstruction.MEASUREMENT:
            qubit = ins.t_qubit[0]
            prob0 = self.get_measurement_probability(qubit, state=0)
            if np.random.rand() < prob0:
                # measured 0
                self.classical_reg.append(0)
                self.single_qubit_unitary(qubit, np.array([[prob0**-0.5, 0], [0, 0]]))
            else:
                # measured 1 
                self.classical_reg.append(1)
                self.single_qubit_unitary(qubit, np.array([[0, 0], [0, (1-prob0)**-0.5]]))
        elif ins.type == GateInstruction.COLLAPSE:
            qubit = ins.t_qubit[0]
            state = ins.params[0]
            prob0 = self.get_measurement_probability(qubit, state=0)
            prob1 = 1 - prob0
            if state == 0:
                mat = np.array([[prob0**-0.5, 0], [0, 0]]) if prob0 > prob1 else np.array([[0, prob1**-0.5], [0, 0]])
            else: # state == 1
                mat = np.array([[0, 0], [0, prob1**-0.5]]) if prob1 > prob0 else np.array([[0, 0], [prob0**-0.5, 0]])
            
            self.single_qubit_unitary(qubit, mat)

    def simulate(self):
        for ins in self.circuit.instructions:
            self.apply_instruction(ins)