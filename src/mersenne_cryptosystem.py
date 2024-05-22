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
    
    # Here we use intp to avoid confusion with python's own int.
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
        '''
        Decodes bit by bit as described in the paper
        '''
        d = ham(self.multiply_bitarray(sk, C))
        dec_bit = ''
        if d <= 2 * (self.h ** 2):
            dec_bit = '0'
        elif d >= self.n - 2 * (self.h ** 2):
            dec_bit = '1'
        else:
            dec_bit = '⊥'
        
        return dec_bit
    
    def KeyGen(self) -> tuple[tuple[bitarray], bitarray]:
        '''
        Generates the public and secret key as is described in the paper.
        '''
        F = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        G = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        R = generate_random_n_bitarray(self.n)
        pk = (R, self.add_bitarray(self.multiply_bitarray(F, R), G))
        sk = F
        return (pk, sk)

    def Enc(self, pk: tuple[bitarray, bitarray], m: bitarray, enc, enc_args: Union[None, dict] = None) -> tuple[bitarray]:
        '''
        Encrypts the message as is described in the paper.

        enc: error correcting function that encodes the message.
        enc_args: the arguments that enc needs.
        '''
        R, T = pk
        A = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        B1 = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        B2 = generate_random_n_hamming_weight_bitarray(self.n, self.h)
        enc_message = enc(m, **enc_args) if enc_args else enc(m)
        # Here we truncate so we can properly bitwise XOR
        return (self.add_bitarray(self.multiply_bitarray(A, R), B1), self.add_bitarray(self.multiply_bitarray(A, T), B2)[-len(enc_message):] ^ enc_message)

    def Dec(self, sk: bitarray, C: tuple[bitarray, bitarray], dec, dec_args: Union[None, dict] = None):
        '''
        Decrypts the message as is described in the paper.

        dec: error correcting function that decodes the message.
        dec_args: the arguments that dec needs.
        '''
        C1, C2 = C
        F = sk
        # Here we truncate so we can properly bitwise XOR
        return dec(self.multiply_bitarray(F, C1)[-len(C2):] ^ C2, **dec_args) if dec_args else dec(self.multiply_bitarray(F, C1)[-len(C2):] ^ C2)
    
    def Encaps(self, pk: tuple[bitarray, bitarray], enc, enc_args: Union[None, dict] = None) -> tuple[tuple[bitarray], bitarray]:
        K = generate_random_n_bitarray(self.lamb)
        '''
        Implements Encaps as is described in the paper.
        '''
        # The expandable hash function is a shuffle where K is the seed used for shuffling.
        # Recommended in https://crypto.stackexchange.com/questions/90156/expandable-hash-function
        random.seed(K.to01())
        R, T = pk
        A = bitarray(random.shuffle([False] * (self.n - self.h) + [True] * self.h))
        B1 = bitarray(random.shuffle([False] * (self.n - self.h) + [True] * self.h))
        B2 = bitarray(random.shuffle([False] * (self.n - self.h) + [True] * self.h))
        enc_K = enc(K, **enc_args) if enc_args else enc(K)
        return ((self.add_bitarray(self.multiply_bitarray(A, R), B1), self.add_bitarray(self.multiply_bitarray(A, T), B2)[-len(enc_K):] ^ enc_K), K)
    
    def Decaps(self, pk: bitarray, sk: bitarray, C: tuple[bitarray, bitarray], enc, dec, enc_args: Union[None, dict] = None, dec_args: Union[None, dict] = None):
        '''
        Implements Decaps as is described in the paper
        '''
        C1, C2 = C
        F = sk
        R, T = pk
        
        K_prime = dec(self.multiply_bitarray(F, C1)[-len(C2):] ^ C2, **dec_args) if dec_args else dec(self.multiply_bitarray(F, C1)[-len(C2):] ^ C2)
        # Again, we use K_prime as the random seed, swapping with K.
        random.seed(K_prime.to01())

        A_prime = bitarray(random.shuffle([False] * (self.n - self.h) + [True] * self.h))
        B1_prime = bitarray(random.shuffle([False] * (self.n - self.h) + [True] * self.h))
        B2_prime = bitarray(random.shuffle([False] * (self.n - self.h) + [True] * self.h))
        enc_K_prime = enc(K_prime, **enc_args) if enc_args else enc(K_prime)
        C_prime = (self.add_bitarray(self.multiply_bitarray(A_prime, R), B1_prime), self.add_bitarray(self.multiply_bitarray(A_prime, T), B2_prime)[-len(enc_K_prime):] ^ enc_K_prime)
        return K_prime if C_prime == C else '⊥'

    
        

    

