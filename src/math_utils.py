import random
from bitarray import bitarray
from bitarray.util import int2ba
from typing import Union

# -------------------- #
#    PRIMALITY TESTS   #
# -------------------- #

def jacobi_symbol(m: int, n: int) -> int:
    if m > n:
        return jacobi_symbol(m % n, n)
    elif m <= 1:
        return m
    elif m % 2 == 0:
        if n % 8 == 1 or n % 8 == 7:
            return jacobi_symbol(m // 2, n)
        else:
            return -jacobi_symbol(m // 2, n)
    else:
        if n % 4 == 3 and m % 4 == 3:
            return -jacobi_symbol(n % m, m)
        else:
            return jacobi_symbol(n % m, m)
        
def solovay_strassen(n: int) -> bool:
    a = random.choice(range(1, n))
    x = jacobi_symbol(a, n)
    if x == 0:
        return True
    y = pow(a, (n - 1) // 2, n)
    if x % n == y:
        return False
    else:
        return True
    
def is_prime(n: int, t: int = 100) -> bool:
    for _ in range(t):
        if solovay_strassen(n):
            return False
    return True

# In this function p is already a Mersenne Number and n is its prime power
def lucas_lehmer(p: int, n: int) -> bool:
    if n == 2:
        return True
    s = 4
    for _ in range(n - 2):
        s = ((s * s) - 2) % p
    return s == 0

def is_mersenne_prime(p: int, n: int) -> bool:
    return is_prime(n) and lucas_lehmer(p, n)


# -------------------- #
#    BIT CONVERTIONS   #
# -------------------- #

def string_to_bits(s: str) -> str:
    return ''.join(format(ord(x), '08b') for x in s)

def bits_to_int(bits: bitarray) -> int:
    result = 0
    for bit in bits:
        result = (result << 1) | bit
    return result

def bits_to_string(bits: bitarray) -> str:
    chars = []
    for i in range(0, len(bits), 8):
        byte = bits[i:i+8]
        chars.append(chr(bits_to_int(byte)))
    return ''.join(chars)

def generate_bit_arrays(n) -> list[bitarray]:
    if n <= 0:
        return [bitarray()]
    else:
        bit_arrays = []
        for bit_array in generate_bit_arrays(n - 1):
            bit_arrays.append(bit_array + [0])
            bit_arrays.append(bit_array + [1])
        return bit_arrays
    
def generate_random_n_bitarray(n: int) -> bitarray:
    return bitarray([random.choice([False, True]) for _ in range(n)])

def generate_random_n_hamming_weight_bitarray(n:int , h: int) -> bitarray:
    string = [False] * (n - h) + [True] * h
    random.shuffle(string)
    return bitarray(string)

# ------------------ #
#    MISCELLANEOUS   #
# ------------------ #

def ham(y: bitarray) -> int:
    return sum(1 for bit in y if bit == True)

# Extended Euclides Algorithm as presented in the textbook
def EED(a: int, b: int) -> tuple[int, int, int]:
    a0 = a
    b0 = b
    t0 = 0
    t = 1
    s0 = 1
    s = 0
    q = a0 // b0 # Entire division
    r = a0 - (q * b0)

    while r > 0:
        temp = t0 - (q * t)
        t0 = t
        t = temp
        temp = s0 - (q * s)
        s0 = s
        s = temp
        a0 = b0
        b0 = r
        q = a0 // b0
        r = a0 - (q * b0)
    r = b0
    return (r, s, t)

def inverse_modulo(element: int, n: int) -> int:
    _, _ ,t = EED(n, element)
    return t % n