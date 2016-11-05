from networkx import Graph


class MoleculeContainer(Graph):
    def __init__(self, meta=None):
        super(MoleculeContainer, self).__init__()
        self.graph['meta'] = meta or {}

    @property
    def meta(self):
        return self.graph['meta']


class ReactionContainer(dict):
    def __init__(self, substrats=None, products=None, meta=None):
        super(ReactionContainer, self).__init__(substrats=substrats or [], products=products or [], meta=meta or {})

    @property
    def substrats(self):
        return self['substrats']

    @property
    def products(self):
        return self['products']

    @property
    def meta(self):
        return self['meta']

    def copy(self):
        return ReactionContainer(substrats=[x.copy() for x in self['substrats']], meta=self['meta'].copy(),
                                 products=[x.copy() for x in self['products']])
