from qiskit.quantum_info import Pauli, SparsePauliOp
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

def xxz_hamiltonian(n, delta=0.5, select=None):
    """Returns the XXZ model Hamiltonian for n qubits and a given delta in Qiskit.
    select = "XX", "YY", "ZZ" or None. None (default) returns XXZ, while the others select
    only the corresponding Pauli terms.
    """
    
    # Initialize lists to store the Pauli strings and coefficients
    pauli_strings = []
    coefficients = []

    for i in range(n):
        # XX term (X_i * X_{i+1})
        x_term = ['I'] * n
        x_term[i] = 'X'
        x_term[(i + 1)%n] = 'X'
        if select == None or select == 'XX':
            pauli_strings.append(''.join(x_term))
            coefficients.append(1.0)
        
        # YY term (Y_i * Y_{i+1})
        y_term = ['I'] * n
        y_term[i] = 'Y'
        y_term[(i + 1)%n] = 'Y'
        if select == None or select == 'YY':
            pauli_strings.append(''.join(y_term))
            coefficients.append(1.0)
        
        # ZZ term (Z_i * Z_{i+1})
        z_term = ['I'] * n
        z_term[i] = 'Z'
        z_term[(i + 1)%n] = 'Z'
        if select == None or select == 'ZZ':
            pauli_strings.append(''.join(z_term))
            coefficients.append(delta)

    # Create the SparsePauliOp object
    paulis = [Pauli(p) for p in pauli_strings]
    hamiltonian = SparsePauliOp(paulis, np.array(coefficients))
    
    return hamiltonian

def binary_code_to_index(key):
    index = 0
    size = len(key)
    for i in range(size):
        index += int(key[i])* 2 ** (size - 1 - i)
    return index

def sample_to_expectation(obs, keys, freq):
    # check observable diagonal
    if (
    np.count_nonzero(
        obs - np.diag(np.diagonal(obs))
    )
    != 0
    ):
        print( "Observable is not diagonal.")
    expval = 0
    for i, k in enumerate(keys):
        index = binary_code_to_index(k)
        expval += obs[index, index] * freq[i]
    return np.real(expval)

def rotate_circuit_XYZ(qc):
    """Generate 
    Args:
        qc: qiskit circuit
    Returns: modified circuits for measuring in computational, 'X', and 'Y' basis.
    """
    nqubits = qc.num_qubits
    # X
    qc_x = deepcopy(qc)
    for i in range(nqubits):
        qc_x.h(i)
    qc_x.measure_all()
    # Y
    qc_y = deepcopy(qc)
    for i in range(nqubits):
        qc_y.sdg(i)
        qc_y.h(i)
    qc_y.measure_all()
    # Z
    qc_z = deepcopy(qc)
    qc_z.measure_all()
    return [qc_x, qc_y, qc_z]


def plot_side_by_side_bars(data_dicts, labels=None, xlabel=None, ylabel=None, title=None):
    """
    Plots side-by-side bar plots for a list of data dictionaries.
    
    Parameters:
    - data_dicts: A list of dictionaries where keys are categories and values are data values.
    
    Each dictionary in the list represents a separate dataset.
    """
    # Combine all keys from the data dictionaries to get a complete list of categories
    all_categories = sorted(set().union(*[data.keys() for data in data_dicts]))
    
    # Create a list of values for each dataset, filling missing categories with 0
    all_values = [[data.get(category, 0) for category in all_categories] for data in data_dicts]
    
    # Number of datasets and categories
    num_datasets = len(data_dicts)
    n = len(all_categories)
    
    # Define the positions and width of the bars
    bar_width = 0.6 / num_datasets  # Dynamic width adjustment based on the number of datasets
    x = np.arange(n)  # the label locations
    
    
    # Plot each dataset
    for i, values in enumerate(all_values):
        if labels is None:
            plt.bar(x + (i - num_datasets / 2) * bar_width + bar_width / 2, values, width=bar_width, label=f'Data Set {i + 1}')
        else:
            plt.bar(x + (i - num_datasets / 2) * bar_width + bar_width / 2, values, width=bar_width, label=labels[i])
    
    # Add labels and title
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(x, all_categories)  # Set the labels for each bar group
    plt.legend()  # Add a legend
    
    # Show the plot
    plt.tight_layout()
    plt.show()
    
def add_partial_swap(qc, qubit1, qubit2, t=np.pi):
    # by default t=pi, the circuit is equivalent to a full swap gate
    # ZZ
    qc.rzz(2*t, qubit1, qubit2)
    qc.barrier()

    # YY 
    qc.s(qubit1) 
    qc.s(qubit2) 
    qc.h(qubit1) 
    qc.h(qubit2) 
    qc.rzz(2*t, qubit1,qubit2)
    qc.h(qubit1) 
    qc.h(qubit2) 
    qc.sdg(qubit1) 
    qc.sdg(qubit2) 

    # XX
    qc.h(qubit1) 
    qc.h(qubit2) 
    qc.barrier()
    qc.rzz(2*t, qubit1,qubit2)
    qc.h(qubit1) 
    qc.h(qubit2)