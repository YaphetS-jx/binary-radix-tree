def print_bits(num: int, num_bits: int) -> None:
    """Print the binary representation of a number up to num_bits."""
    bits = bin(num)[2:].zfill(num_bits)
    print(bits)

def interleave_bits(x: int, y: int, z: int, num_bits_per_dim: int) -> int:
    """
    Interleave bits from x, y, and z coordinates into a single Morton code.
    
    Args:
        x, y, z: Integer coordinates
        num_bits: Number of bits to interleave from each coordinate
    
    Returns:
        Interleaved Morton code
    """    
    
    result = 0
    for i in range(num_bits_per_dim):
        result |= ((x & (1 << i)) << (2 * i)) | \
                 ((y & (1 << i)) << (2 * i + 1)) | \
                 ((z & (1 << i)) << (2 * i + 2))
    return result

def deinterleave_bits(morton_code: int, num_bits_per_dim: int) -> tuple[int, int, int]:
    """
    Deinterleave bits from a Morton code into x, y, z coordinates.
    
    Args:
        morton_code: The Morton code to deinterleave
        num_bits_per_dim: Number of bits per dimension
    
    Returns:
        Tuple of (x, y, z) coordinates
    """    
    
    x = y = z = 0
    for i in range(num_bits_per_dim):
        x |= (morton_code & (1 << (3 * i))) >> (2 * i)
        y |= (morton_code & (1 << (3 * i + 1))) >> (2 * i + 1)
        z |= (morton_code & (1 << (3 * i + 2))) >> (2 * i + 2)
    return x, y, z

def fractional_to_binary(num: float, num_bits: int) -> int:
    """
    Convert a fractional number between 0 and 1 to binary representation.
    
    Args:
        num: Number between 0 and 1
        num_bits: Number of bits to use
    
    Returns:
        Binary representation as integer
    """
    assert 0 <= num < 1.0
    
    result = 0
    n = 0
    
    while num > 0 and n < num_bits:
        num *= 2
        if num >= 1.0:
            result |= 1 << (num_bits - 1 - n)
            num -= 1.0
        n += 1
    return result

def coordinate_to_morton_code(num: float, num_bits: int) -> int:
    """
    Convert a coordinate (0 to 1) to a Morton code.
    
    Args:
        num: Coordinate value between 0 and 1
        num_bits: Number of bits to use
    
    Returns:
        Morton code for the coordinate
    """
    assert 0 <= num < 1.0
    return fractional_to_binary(num, num_bits - 1)

def binary_to_fractional(binary: int, num_bits: int) -> float:
    """
    Convert a binary number to a fractional value between 0 and 1.
    
    Args:
        binary: Binary number
        num_bits: Number of bits to consider
    
    Returns:
        Fractional value between 0 and 1
    """
    assert num_bits <= 64
    
    result = 0.0
    for i in range(num_bits):
        result += (binary & 1) * 1.0 / (1 << (num_bits - i))
        binary >>= 1
    return result

def morton_code_to_coordinate(morton_code: int, num_bits: int) -> float:
    """
    Convert a Morton code back to a coordinate value.
    
    Args:
        morton_code: The Morton code
        num_bits: Number of bits used in the Morton code
    
    Returns:
        Coordinate value between 0 and 1
    """
    assert num_bits >= 1 and 3 * num_bits <= 64
    assert (morton_code & (1 << (num_bits - 1))) == 0
    
    return binary_to_fractional(morton_code, num_bits - 1)


def coordinate_to_morton_code_3d(x: float, y: float, z: float, num_bits_per_dim: int) -> int:
    """
    Convert a coordinate (0 to 1) to a Morton code.
    
    Args:
        x, y, z: Coordinate values between 0 and 1
        num_bits_per_dim: Number of bits per dimension
    
    Returns:
        Morton code for the coordinate
    """
    x_mc = coordinate_to_morton_code(x, num_bits_per_dim)
    y_mc = coordinate_to_morton_code(y, num_bits_per_dim)
    z_mc = coordinate_to_morton_code(z, num_bits_per_dim)
    
    return interleave_bits(x_mc, y_mc, z_mc, num_bits_per_dim)


