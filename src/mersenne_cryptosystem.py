import sys
import os

notebook_dir = os.getcwd()
parent_dir = os.path.dirname(notebook_dir)
src_dir = os.path.join(parent_dir, 'src')

sys.path.append(src_dir)

from typing import Literal
from math_utils import *
from reed_muller import *

class MersenneCryptosystem():
    """
    A class representing a Mersenne Cryptosystem.

    Args:
        lamb (int): Security Parameter.
        n (int): Exponent of Mersenne Prime.
        h (int): Low Hamming Weight.

    Fields:
        p (int): The Mersenne prime
    """
    def __init__(self, lamb: int, n: int, h: int) -> None:
        self.lamb = lamb
        self.p = (2 ** n) - 1
        self.n = n
        self.h = h

    def seq(self, x: int) -> bitarray:
        '''
        Converts integers into bit arrays representing their representative in Zp.
        '''
        x = x % self.p
        return int2ba(x, self.n)
    
    def intp(self, y: bitarray) -> int:
        '''
        Converts bitarrays into the representative of their integer in Zp.
        '''
        return bits_to_int(y) % self.p
    
    def add_bitarray(self, y1: bitarray, y2: bitarray) -> bitarray:
        '''
        Addition is performed by first calculating the integer representation in Zp of the bitarrays, adding them up and then calculating their bitarray representation in Zp.
        '''
        added_int = self.intp(y1) + self.intp(y2)
        return self.seq(added_int)
    
    def multiply_bitarray(self, y1: Union[bitarray, int], y2: bitarray) -> bitarray:
        '''
        Multiplication is performed by first calculating the integer representation in Zp of the bitarrays, adding them up and then calculating their bitarray representation in Zp.
        '''
        if y1 == -1:
            multiplied_int = -1 * self.intp(y2)
        elif y1 == 1:
            multiplied_int = self.intp(y2)
        else:
            multiplied_int = self.intp(y1) * self.intp(y2)
        return self.seq(multiplied_int)
    
    def divide_bitarray(self, y1: bitarray, y2: bitarray) -> bitarray:
        '''
        Division is performed by first calculating the inverse in Zp of the integer representation in Zp of y2 and then multiplying their integer representation in Zp in mod p.
        '''
        inv_y2 = inverse_modulo(self.intp(y2), self.p)
        return (self.intp(y1) * inv_y2) % self.p

    def bit_by_bit_key_gen(self) -> tuple[bitarray]:
        '''
        Generates the keys for bit by bit encryption as described in the paper.
        '''
        F = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        G = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        pk = self.seq(self.divide_bitarray(F, G))
        sk = G
        return (pk, sk)

    def bit_by_bit_enc(self, pk: bitarray, b: int) -> bitarray:
        '''
        Encodes a bit as described in the paper.
        '''
        A = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        B = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        if b == 0:
            b = 1
        else:
            b = -1
        return self.multiply_bitarray(b, self.add_bitarray(self.multiply_bitarray(A, pk), B))
        

    def bit_by_bit_dec(self, sk: bitarray, C: bitarray) -> Literal['0', '1', '⊥']:
        d = ham(self.multiply_bitarray(sk, C))
        dec_bit = ''
        if d <= 2 * (self.h ** 2):
            dec_bit = '0'
        elif d >= self.n - 2 * (self.h ** 2):
            dec_bit = '1'
        else:
            dec_bit = '⊥'
        
        return dec_bit
        

    

