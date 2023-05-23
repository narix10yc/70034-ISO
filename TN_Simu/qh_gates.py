import numpy as np 

### Single-qubit Non-parametric ###
def _get_X():
    return np.array([[0, 1], [1, 0]])

def _get_Y():
    return np.array([[0, -1j], [1j, 0]])

def _get_Z():
    return np.array([[1, 0], [0, -1]])

def _get_H():
    return 2**-0.5 * np.array([[1, 1], [1, -1]])

def _get_S():
    return np.array([[1, 0], [0, 1j]])

def _get_T():
    a = 2 ** -0.5
    return np.array([[1, 0], [0, a+a*1j]])


### Single-qubit Parametric ###
def _get_P(theta):
    t = np.exp(theta*1j)
    return np.array([[1, 0], [0, t]])

def _get_RX(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -1j*s], [-1j*s, c]])

def _get_RY(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -s], [s, c]])

def _get_RZ(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c-s*1j, 0], [0, c+s*1j]])

def _get_U(theta, phi, lam):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -np.exp(lam*1j) * s], [np.exp(phi*1j) * s, np.exp((phi+lam)*1j) * c]])


### Two-qubit Non-parametric ###
def _get_CX():
    return np.array([[1,0,0,0],  [0,0,0,1],  [0,0,1,0],  [0,1,0,0]]).reshape(2,2,2,2)

def _get_CZ():
    return np.array([[1,0,0,0],  [0,1,0,0],  [0,0,1,0],  [0,0,0,-1]]).reshape(2,2,2,2)


### Two-qubit Parametric ###
def _get_RXX(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([    c,     0,     0, -1j*s, 
                         0,     c, -1j*s,     0,
                         0, -1j*s,     c,     0,
                     -1j*s,     0,     0,     c]).reshape(2,2,2,2)

def _get_RYY(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([    c,     0,     0,  1j*s, 
                         0,     c, -1j*s,     0,
                         0, -1j*s,     c,     0,
                      1j*s,     0,     0,     c]).reshape(2,2,2,2)

def _get_RZZ(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([c-1j*s,0,0,0, 
                     0,c+1j*s,0,0, 
                     0,0,c+1j*s,0, 
                     0,0,0,c-1j*s]).reshape(2,2,2,2)



### Public Methods ###
gate_dict = {
    # single-qubit non-para
    "X": _get_X,
    "Y": _get_Y,
    "Z": _get_Z,
    "H": _get_H, 
    "S": _get_S,
    "T": _get_T,
    # single-qubit para
    "P": _get_P,
    "RX": _get_RX,
    "RY": _get_RY,
    "RZ": _get_RZ,
    "U": _get_U,
    # two-qubit non-para
    "CX": _get_CX,
    "CZ": _get_CZ,
    # two-qubit para
    "RXX": _get_RXX,
    "RYY": _get_RYY,
    "RZZ": _get_RZZ,

}



def get_matrix(gate: str, *params):
    f = gate_dict[gate.upper()]
    return f(*params)
