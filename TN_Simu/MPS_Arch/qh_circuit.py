import numpy as np 

class GateInstruction():
    UNITARY_GATE = 0
    MEASUREMENT = 1
    COLLAPSE = 2
    def __init__(self, gate, t_qubit, c_qubit=[], *params):
        self.gate = gate 
        self.t_qubit = t_qubit
        self.c_qubit = c_qubit
        self.params = params
        self.nqubits = len(self.t_qubit) + len(self.c_qubit)
        self.type = GateInstruction.UNITARY_GATE
    
    def __str__(self):
        s = f"{self.gate}"
        if len(self.params) > 0:
            s += "(" + ' '.join(['{:.3f}'.format(param) for param in self.params]) + ")"
        s += f" on {self.t_qubit}"
        if len(self.c_qubit) > 0:
            s += f" controled on {self.c_qubit}"
        return s 


class QHCircuit():
    def __init__(self, nqubits):
        self.nqubits: int = nqubits
        self.instructions: list[GateInstruction] = []
    
    def __str__(self):
        return '\n'.join([ins.__str__() for ins in self.instructions])

    def _add_instruction(self, ins: GateInstruction):
        self.instructions.append(ins)

    def x(self, qubit):
        self._add_instruction(GateInstruction("X", [qubit]))

    def y(self, qubit):
        self._add_instruction(GateInstruction("Y", [qubit]))

    def z(self, qubit):
        self._add_instruction(GateInstruction("Z", [qubit]))

    def h(self, qubit):
        self._add_instruction(GateInstruction("H", [qubit]))

    def s(self, qubit):
        self._add_instruction(GateInstruction("S", [qubit]))

    def t(self, qubit):
        self._add_instruction(GateInstruction("T", [qubit]))

    def cx(self, control, target):
        self._add_instruction(GateInstruction("CX", [target], [control]))

    def cz(self, control, target):
        self._add_instruction(GateInstruction("CZ", [target], [control]))

    def cnot(self, control, target):
        self.cx(control, target)
    

    def p(self, qubit, theta):
        self._add_instruction(GateInstruction("P", [qubit], [], theta))

    def rx(self, qubit, theta):
        self._add_instruction(GateInstruction("RX", [qubit], [], theta))

    def ry(self, qubit, theta):
        self._add_instruction(GateInstruction("RY", [qubit], [], theta))

    def rz(self, qubit, theta):
        self._add_instruction(GateInstruction("RZ", [qubit], [], theta))
    
    def u(self, qubit, theta, phi, lam):
        self._add_instruction(GateInstruction("U", [qubit], [], theta, phi, lam))


    def rxx(self, qubit1, qubit2, theta):
        self._add_instruction(GateInstruction("RXX", [qubit1, qubit2], [], theta))

    def ryy(self, qubit1, qubit2, theta):
        self._add_instruction(GateInstruction("RYY", [qubit1, qubit2], [], theta))

    def rzz(self, qubit1, qubit2, theta):
        self._add_instruction(GateInstruction("RZZ", [qubit1, qubit2], [], theta))


    def measure(self, qubit):
        ins = GateInstruction("MEASURE", [qubit])
        ins.type = GateInstruction.MEASUREMENT
        self._add_instruction(ins)
    
    def measure_all(self):
        for qubit in range(self.nqubits):
            self.measure(qubit)
    
    def collapse(self, qubit, state: int):
        ins = GateInstruction("COLLAPSE", [qubit], [], int(state))
        ins.type = GateInstruction.COLLAPSE
        self._add_instruction(ins)












