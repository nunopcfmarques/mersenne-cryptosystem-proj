from bitarray import bitarray

def enc_to_repetition_code(bits: bitarray, p: int) -> bitarray:
    repeated_bits = bitarray()
    for bit in bits:
        repeated_bits.extend([bit] * p)
    return repeated_bits

def dec_repetition_code_to_bits(enc_bits: bitarray, p: int) -> bitarray:
    decoded_bits = bitarray()
    for i in range(0, len(enc_bits), p):
        chunk = enc_bits[i:i+p]
        if chunk.count(True) > p / 2:
            decoded_bits.append(True)
        else:
            decoded_bits.append(False)
    return decoded_bits