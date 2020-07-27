

def augmented(molecule, deep):
    bonds = molecule._bonds

    if deep < 1:
        raise ValueError('Deep should be >= 1')
    stack = []
    groups = set()
    response = []
    deep -= 1
    for a, n in bonds.items():
        n = [x for x in n]
        aug = (a, *n)
        if aug not in groups:
            groups.add(aug)
            response.append(molecule.substructure(aug, as_query=True))
            stack.append((list(aug), [*n], deep))

    while stack:
        augs, nei, st = stack.pop(0)
        st -= 1
        level = []
        for a in nei:
            n = list(bonds[a])
            level.extend(n)
        augs.extend(level)
        aug = tuple(set(augs))
        if aug not in groups:
            groups.add(aug)
            response.append(molecule.substructure(aug, as_query=True))
            if st:
                stack.append((augs, level, st))
    return response
