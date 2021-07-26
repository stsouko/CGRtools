

def unpack(bytes data):
    cdef short isotope_shift
    cdef unsigned char a, b, c, d
    cdef unsigned short na, nct, i, n, shift = 3
    cdef unsigned long bb, nb = 0

    cdef unsigned short[4095] atom, neighbors, hydrogens, orders, mapping, isotopes, cis_trans_1, cis_trans_2
    cdef unsigned short[8190] connections
    cdef short[4095] charges
    cdef bint[4095] radicals, is_tet, is_all, tet_sign, all_sign, ct_sign
    cdef float[4095] x, y

    cdef dict py_charges, py_radicals, py_hydrogens, py_plane
    cdef dict py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo
    cdef tuple py_xy
    cdef list py_mapping, py_atoms, py_isotopes, py_neighbors, py_connections, py_orders

    a, b, c = data[:3]
    na = a << 4| (b & 0xf0) >> 4
    nct = (b & 0x0f) << 8 | c

    for i in range(na):
        a, b = data[shift: shift + 2]
        mapping[i] = a << 4 | (b & 0xf0) >> 4
        neighbors[i] = b & 0x0f
        nb += b & 0x0f

        a, b = data[shift + 2: shift + 4]
        if a & 0x80:
            is_tet[i] = True
            tet_sign[i] = a & 0x40
        else:
            is_tet[i] = False
        if a & 0x20:
            is_all[i] = True
            all_sign[i] = a & 0x10
        else:
            is_all[i] = False

        atom[i] = b & 0x7f
        isotope_shift = (a & 0x0f) << 1 | b >> 7
        if isotope_shift:
            isotopes[i] = common_isotopes[b & 0x7f] + isotope_shift
        else:
            isotopes[i] = 0

        a, b = data[shift + 4: shift + 6]
        x[i] = a << 8 | b
        a, b = data[shift + 6: shift + 8]
        y[i] = a << 8 | b

        a = data[shift + 8]
        hydrogens[i] = a >> 5
        charges[i] = ((a >> 1) & 0x0f) - 4
        radicals[i] = a & 0x01

        shift += 9

    nb //= 2
    for i in range(nb):
        a, b, c = data[shift: shift + 3]
        connections[i * 2] = a << 4| (b & 0xf0) >> 4
        connections[i * 2 + 1] = (b & 0x0f) << 8 | c
        shift += 3

    if nb % 5:
        bb = nb // 5 + 1
    else:
        bb = nb // 5

    for i in range(bb):
        a, b = data[shift: shift + 2]
        orders[i * 5] = (a >> 4) + 1
        orders[i * 5 + 1] = ((a >> 1) & 0x07) + 1
        orders[i * 5 + 2] = ((a & 0x01) << 2 | b >> 6) + 1
        orders[i * 5 + 3] = ((b >> 3) & 0x07) + 1
        orders[i * 5 + 4] = (b & 0x07) + 1
        shift += 2

    for i in range(nct):
        a, b, c, d = data[shift: shift + 4]
        cis_trans_1[i] = a << 4 | (b & 0xf0) >> 4
        cis_trans_2[i] = (b & 0x0f) << 8 | c
        ct_sign[i] = d
        shift += 4

    py_mapping = []
    py_atoms = []
    py_isotopes = []
    py_neighbors = []
    py_connections = []
    py_orders = []

    py_charges = {}
    py_radicals = {}
    py_hydrogens = {}
    py_plane = {}
    py_atoms_stereo = {}
    py_allenes_stereo = {}
    py_cis_trans_stereo = {}
    for i in range(na):
        n = mapping[i]
        py_mapping.append(n)
        py_atoms.append(atom[i])
        py_isotopes.append(isotopes[i])
        py_neighbors.append(neighbors[i])

        py_charges[n] = charges[i]
        py_radicals[n] = radicals[i]
        py_hydrogens[n] = hydrogens[i]

        py_xy = (x[i], y[i])
        py_plane[n] = py_xy

        if is_tet[i]:
            py_atoms_stereo[n] = tet_sign[i]
        if is_all[i]:
            py_allenes_stereo[n] = all_sign[i]

    for i in range(nb):
        py_orders.append(orders[i])
    for i in range(nb * 2):
        py_connections.append(connections[i])

    for i in range(nct):
        py_xy = (cis_trans_1[i], cis_trans_2[i])
        py_cis_trans_stereo[py_xy] = ct_sign[i]

    return (py_mapping, py_atoms, py_isotopes, py_neighbors, py_connections, py_orders,
            py_charges, py_radicals, py_hydrogens, py_plane,
            py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo)


cdef short[119] common_isotopes
common_isotopes[:] = [0, -15, -12, -9, -7, -5, -4, -2, 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 24, 23, 24, 29,
                      32, 35, 36, 39, 40, 43, 43, 48, 49, 54, 57, 59, 63, 64, 68, 69, 72, 73, 75, 77,
                      80, 82, 85, 87, 90, 92, 96, 99, 103, 106, 112, 111, 115, 117, 121, 123, 124, 125,
                      128, 129, 134, 136, 141, 143, 147, 149, 151, 153, 157, 159, 162, 165, 168, 170,
                      174, 176, 179, 181, 185, 188, 191, 193, 193, 194, 206, 207, 210, 211, 216, 215,
                      222, 221, 228, 227, 231, 231, 235, 236, 241, 242, 243, 244, 245, 254, 253, 254,
                      254, 262, 265, 265, 269, 262, 273, 273, 277, 281, 278]


from timeit import timeit
