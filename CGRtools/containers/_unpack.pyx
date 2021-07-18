

def unpack(const unsigned char[:] data):
    cdef unsigned char a, b, c, d
    cdef unsigned short na, nct

    a, b, c = data[:3]
    na = a << 4| (b & 0xf0) >> 4
    nct = (b & 0x0f) << 4 | c

    return na, nct
