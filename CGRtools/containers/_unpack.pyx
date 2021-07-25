

def unpack(const unsigned char[:] data):
    cdef short isotope_shift
    cdef unsigned char a, b, c, d
    cdef unsigned short na, nct, tet, all, i, shift = 3

    cdef unsigned char[4095] atom, neighbors, hydrogens
    cdef unsigned short[4095] mapping, isotope
    cdef char[4095] charge
    cdef bint[4095] radical, is_tet, is_all, tet_sign, all_sign
    cdef float[4095] x, y

    a, b, c = data[:3]
    na = a << 4| (b & 0xf0) >> 4
    nct = (b & 0x0f) << 8 | c

    for i in range(na):
        a, b = data[shift: shift + 2]
        mapping[i] = a << 4 | (b & 0xf0) >> 4
        neighbors[i] = b & 0x0f

        a, b = data[shift + 2: shift + 4]
        if a & 0x80:
            is_tet[i] = True
            tet_sign[i] = a & 0x40
        if a & 0x20:
            is_all[i] = True
            all_sign[i] = a & 0x10

        atom[i] = b & 0x7f
        isotope_shift = (a & 0x0f) << 1 | b >> 7
        if isotope_shift:
            isotope[i] = common_isotopes[b & 0x7f] + isotope_shift
        else:
            isotope[i] = 0

        a, b = data[shift + 4: shift + 6]
        x[i] = a << 8 | b
        a, b = data[shift + 6: shift + 8]
        y[i] = a << 8 | b

        a = data[shift + 8]
        hydrogens[i] = a >> 5
        charge[i] = ((a >> 1) & 0x0f) - 4
        radical[i] = a & 0x01

        shift += 9

    return


cdef short[119] common_isotopes
common_isotopes[:] = [0, -15, -12, -9, -7, -5, -4, -2, 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 24, 23, 24, 29,
                                32, 35, 36, 39, 40, 43, 43, 48, 49, 54, 57, 59, 63, 64, 68, 69, 72, 73, 75, 77,
                                80, 82, 85, 87, 90, 92, 96, 99, 103, 106, 112, 111, 115, 117, 121, 123, 124, 125,
                                128, 129, 134, 136, 141, 143, 147, 149, 151, 153, 157, 159, 162, 165, 168, 170,
                                174, 176, 179, 181, 185, 188, 191, 193, 193, 194, 206, 207, 210, 211, 216, 215,
                                222, 221, 228, 227, 231, 231, 235, 236, 241, 242, 243, 244, 245, 254, 253, 254,
                                254, 262, 265, 265, 269, 262, 273, 273, 277, 281, 278]
