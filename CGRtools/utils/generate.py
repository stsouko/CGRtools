

def augmented(molecule, deep):
    bonds = molecule._bonds

    if deep < 1:
        raise ValueError('Deep should be >= 1')

    response = []
    groups = set()
    stack = [([a], list(n)) for a, n in bonds.items()]
    while stack:
        aug, nei = stack.pop(0)
        for x in nei:
            augx = (*aug, x)
            if augx not in groups:
                groups.add(augx)
                response.append(molecule.substructure(augx, as_query=True))
                nt = nei.copy()
                nt.remove(x)
                nt.extend(list(bonds[x]))
                if len(augx) < deep:
                    stack.append((augx, nt))

    return response
