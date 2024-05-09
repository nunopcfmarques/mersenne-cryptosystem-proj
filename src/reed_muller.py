from math_utils import *
import itertools
from typing import Iterable

def dot_product(x: bitarray, y: bitarray) -> int:   
    result = sum(a & b for a, b in zip(x, y))
    return result

def matrix_product(x: Union[list[bitarray], bitarray], y: Union[list[bitarray], bitarray]) -> list[bitarray]:
    '''
    Bitwise matrix product of bitarrays or lists of bitarrays. The bitwise matrix product is the usual matrix product where the dot product is calculated modulo 2.
    '''
    result = []
    if isinstance(x, bitarray):
        x = [x]
    if isinstance(y, bitarray):
        y = [y]
    for i in range(len(x)):
        row_result = bitarray()
        for j in range(len(y[0])):
            column_bits = bitarray()
            for row in y:
                column_bits.append(row[j])
            dot_prod = dot_product(x[i], column_bits) % 2
            row_result.append(dot_prod)
        result.append(row_result)
    return result if len(result) > 1 else result[0]

def add_bitarrays(x: bitarray, y: bitarray):
    '''
    Bit arrays are added modulo 2.
    '''
    return bitarray([((a + b) % 2) for a, b in zip(x, y)])

class ReedMullerCode():
    """
    A class representing a ReedMullerCode.

    Args:
        r (int): Degree at most r.
        m (int): amount of variables.

    Fields:
        variables (tuple[str]): Tuple containing the variables from x0 to xm.
        generator_matrix (list[bitarray]): The generator matrix as described in 'Reed-Muller Codes'.
        k (int): Size of messages that can be encoded.
        n (int): Size of encoded message.
        max_erros (int): Maximum errors the code can correct.
    """

    def __init__(self, r: int, m: int) -> None:
        self.r = r
        self.m = m
        self.variables = [f"x{i}" for i in range(self.m)]
        self.generator_matrix, self.degree_to_rows = self.construct_matrix()
        self.k = len(self.generator_matrix) # length of message it can encode
        self.n = 2 ** m # length of codeword
        self.max_errors = 2 ** (self.m - self.r - 1) - 1

    def phi(self, var: str) -> bitarray:
        '''
        Phi function defined in the report titled 'Reed-Muller Codes', more specifically in section 2.
        '''
        index = int(var[-1]) #index will always be last character
        sequence_length = 2 ** self.m
        repeat_times = 2 ** index
        segment_length = sequence_length // (2 ** (index + 1))
        pattern = bitarray([1] * segment_length + [0] * segment_length)

        result = pattern * repeat_times
        return ~result if var[0] == "-" else result
        
    def row_phi(self, subset: tuple[str]) -> bitarray:
        '''
        The phi of a row is the bitwise product of the phi's of its elements
        '''
        row_phi = bitarray([1] * (2 ** self.m))
        for var in subset:
            row_phi &= self.phi(var)
        return row_phi
        
    def subsets(self, i: int) -> Iterable[tuple[str]]:
        '''
        Gives all the possible subsets of size i for the defined variables
        '''
        return itertools.combinations(self.variables, i)
    
    def characteristic_vectors(self, subset: tuple[str]) -> list[tuple[str]]:
        '''
        Calculates the characteristic vectors according to the algorithm defined in the report titled 'Reed Muller Codes', in section 3.1.
        '''
        if subset == ("dummy"):
            subset = ()
        variables_not_in_subset = set(self.variables) - set(subset)
        J = set()

        for variable in variables_not_in_subset:
            J.add("-" + variable)
            J.add(variable)
        degree = self.m - len(subset)
        combos = itertools.combinations(J, degree)
        result = []
        
        # Removes combinations that have both the variable and its complement
        for combo in combos:
            indexes = {var[-1] for var in combo}
            if len(indexes) == len(combo):
                result.append(combo)
        
        return result
    
    def construct_matrix(self) -> list[bitarray]:
        '''
        Implements the algorithm to build a generator matrix defined in the report titled 'Reed Muller Codes', more specifically in section 3.1.
        '''
        matrix = [bitarray([1] * (2 ** self.m))]
        # dummy is used because otherwise the algorithm won't execute 'for row in self.degree_to_rows[degree][::-1]':
        degree_to_rows = {0 : ("dummy",)}
        curr_index = 0
        for i in range(1, self.r + 1):
            subsets = self.subsets(i)
            for subset in subsets:
                curr_index += 1
                matrix.append(self.row_phi(subset))
                degree = len(subset)
                if degree in degree_to_rows:
                    degree_to_rows[degree].append(subset)
                else:
                    degree_to_rows[degree] = [subset]
        return matrix, degree_to_rows
    
    def encode(self, message: bitarray) -> bitarray:
        '''
        The encoding is the matrix product between the message and the generator matrix
        '''
        return matrix_product(message, self.generator_matrix)
    
    def decode(self, message: bitarray) -> bitarray:
        '''
        Implements the majority logic algorithm defined in the report titled 'Reed Muller Codes', more specifically in section 3.2.
        '''
        decoded_message = bitarray(bitarray([0]) * (len(self.generator_matrix)))
        degree = self.r
        pos = len(self.generator_matrix) - 1

        while degree != -1:
            for row in self.degree_to_rows[degree][::-1]:
                count_0 = 0
                count_1 = 0
                for char_vector in self.characteristic_vectors(row):
                    res = dot_product(message, self.row_phi(char_vector)) % 2
                    if res == 1:
                        count_1 += 1
                    else:
                        count_0 += 1
                if count_0 != count_1:
                    decoded_message[pos] = 0 if count_0 > count_1 else 1
                    pos -= 1
            s = matrix_product(decoded_message[pos + 1: pos + len(self.degree_to_rows[degree]) + 1], self.generator_matrix[pos + 1: pos + len(self.degree_to_rows[degree]) + 1])
            message = add_bitarrays(s, message)
            degree -= 1

        decoded_message += bitarray([0]) * (len(self.generator_matrix) - len(decoded_message))
        return decoded_message