# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from CachedMethods import FrozenDict
from .element import Element
from .groups import GroupIII
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Sc(Element, PeriodIV, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 21

    @property
    def isotopes_distribution(self):
        return FrozenDict({45: 1.0, 44: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({45: 44.95591, 44: 43.959403})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()), (-3, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F')))

    @property
    def atomic_radius(self):
        return 1.84


class Y(Element, PeriodV, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 39

    @property
    def isotopes_distribution(self):
        return FrozenDict({89: 1.0, 86: 0., 90: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({89: 88.905848, 86: 85.914886, 90: 89.907152})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.12


class La(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 57

    @property
    def isotopes_distribution(self):
        return FrozenDict({138: 0.0009, 139: 0.9991})

    @property
    def isotopes_masses(self):
        return FrozenDict({138: 137.907107, 139: 138.906348})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.12  # unknown, taken radius of previous element in group


class Ce(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 58

    @property
    def isotopes_distribution(self):
        return FrozenDict({136: 0.00185, 138: 0.00251, 140: 0.8845, 142: 0.11114})

    @property
    def isotopes_masses(self):
        return FrozenDict({136: 135.90714, 138: 137.905986, 140: 139.905434, 142: 141.90924})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.12  # unknown, taken radius of previous element in group


class Pr(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 59

    @property
    def isotopes_distribution(self):
        return FrozenDict({141: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({141: 140.907648})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.47


class Nd(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 60

    @property
    def isotopes_distribution(self):
        return FrozenDict({142: 0.272, 143: 0.122, 144: 0.238, 145: 0.083, 146: 0.172, 148: 0.057, 150: 0.056})

    @property
    def isotopes_masses(self):
        return FrozenDict({142: 141.907719, 143: 142.90981, 144: 143.910083, 145: 144.912569, 146: 145.913112,
                           148: 147.916889, 150: 149.920887})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'),)),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))))

    @property
    def atomic_radius(self):
        return 2.06


class Pm(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 61

    @property
    def isotopes_distribution(self):
        return FrozenDict({145: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({145: 144.912749})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.05


class Sm(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 62

    @property
    def isotopes_distribution(self):
        return FrozenDict({144: 0.0307, 147: 0.1499, 148: 0.1124, 149: 0.1382, 150: 0.0738, 152: 0.2675, 154: 0.2275,
                           145: 0., 153: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({144: 143.911995, 147: 146.914893, 148: 147.914818, 149: 148.91718, 150: 149.917271,
                           152: 151.919728, 154: 153.922205, 145: 144.913410, 153: 152.922097})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.38


class Eu(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 63

    @property
    def isotopes_distribution(self):
        return FrozenDict({151: 0.4781, 153: 0.5219, 152: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({151: 150.919846, 153: 152.921226, 152: 151.921744})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.31


class Gd(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 64

    @property
    def isotopes_distribution(self):
        return FrozenDict({152: 0.002, 154: 0.0218, 155: 0.148, 156: 0.2047, 157: 0.1565, 158: 0.2484, 160: 0.2186,
                           153: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({152: 151.919788, 154: 153.920862, 155: 154.922619, 156: 155.92212, 157: 156.923957,
                           158: 157.924101, 160: 159.927051, 153: 152.921750})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.33


class Tb(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 65

    @property
    def isotopes_distribution(self):
        return FrozenDict({159: 1.0, 160: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({159: 158.925343, 160: 159.927168})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.25


class Dy(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 66

    @property
    def isotopes_distribution(self):
        return FrozenDict({156: 0.0006, 158: 0.001, 160: 0.0234, 161: 0.1891, 162: 0.2551, 163: 0.249, 164: 0.2818})

    @property
    def isotopes_masses(self):
        return FrozenDict({156: 155.924278, 158: 157.924405, 160: 159.925194, 161: 160.92693, 162: 161.926795,
                           163: 162.928728, 164: 163.929171})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.28


class Ho(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 67

    @property
    def isotopes_distribution(self):
        return FrozenDict({165: 1.0, 166: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({165: 164.930319, 166: 165.932284})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.26


class Er(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 68

    @property
    def isotopes_distribution(self):
        return FrozenDict({162: 0.0014, 164: 0.0161, 166: 0.3361, 167: 0.2293, 168: 0.2678, 170: 0.1493})

    @property
    def isotopes_masses(self):
        return FrozenDict({162: 161.928775, 164: 163.929197, 166: 165.93029, 167: 166.932045, 168: 167.932368,
                           170: 169.93546})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.26


class Tm(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 69

    @property
    def isotopes_distribution(self):
        return FrozenDict({169: 1.0, 170: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({169: 168.934211, 170: 169.935801})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.22


class Yb(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 70

    @property
    def isotopes_distribution(self):
        return FrozenDict({168: 0.0013, 170: 0.0304, 171: 0.1428, 172: 0.2183, 173: 0.1613, 174: 0.3183, 176: 0.1276,
                           169: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({168: 167.933894, 170: 169.934759, 171: 170.936322, 172: 171.936378, 173: 172.938207,
                           174: 173.938858, 176: 175.942568, 169: 168.935190})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.22


class Lu(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 71

    @property
    def isotopes_distribution(self):
        return FrozenDict({175: 0.9741, 176: 0.0259, 177: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({175: 174.940768, 176: 175.942682, 177: 176.943758})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17


class Ac(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 89

    @property
    def isotopes_distribution(self):
        return FrozenDict({227: 1.0, 225: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({227: 227.027752, 225: 225.023230})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Th(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 90

    @property
    def isotopes_distribution(self):
        return FrozenDict({232: 1.0, 227: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({232: 232.03805, 227: 227.027704})

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return ((4, False, 0, ()),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Pa(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 91

    @property
    def isotopes_distribution(self):
        return FrozenDict({231: 1.0, 233: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({231: 231.035879, 233: 233.040247})

    @property
    def _common_valences(self):
        return 0, 4, 5

    @property
    def _valences_exceptions(self):
        return ((4, False, 0, ()),
                (0, False, 0, ((1, 'H'), (1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class U(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 92

    @property
    def isotopes_distribution(self):
        return FrozenDict({234: 5.5e-05, 235: 0.0072, 238: 0.992745})

    @property
    def isotopes_masses(self):
        return FrozenDict({234: 234.040946, 235: 235.043923, 238: 238.050783})

    @property
    def _common_valences(self):
        return 0, 3, 4, 5, 6

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (4, False, 0, ()),
                (2, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Np(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 93

    @property
    def isotopes_distribution(self):
        return FrozenDict({237: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({237: 237.048173})

    @property
    def _common_valences(self):
        return 0, 2, 3, 4, 5, 6, 7

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (4, False, 0, ()),
                (1, False, 0, ((2, 'O'), (2, 'O'))),
                (2, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Pu(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 94

    @property
    def isotopes_distribution(self):
        return FrozenDict({239: 1.0, 242: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({239: 239.052163, 242: 242.058743})

    @property
    def _common_valences(self):
        return 0, 3, 4, 5, 6

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (4, False, 0, ()),
                (1, False, 0, ((2, 'O'), (2, 'O'))),
                (2, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'Se'), )),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((2, 'Te'),)),
                (0, False, 0, ((2, 'O'),)),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Am(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 95

    @property
    def isotopes_distribution(self):
        return FrozenDict({241: 1.0, 243: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({241: 241.056829, 243: 243.061380})

    @property
    def _common_valences(self):
        return 0, 2, 3, 4

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Cm(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 96

    @property
    def isotopes_distribution(self):
        return FrozenDict({244: 1.0, 243: 0., 248: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({244: 244.062753, 243: 243.061389, 248: 248.072349})

    @property
    def _common_valences(self):
        return 0, 3, 4

    @property
    def _valences_exceptions(self):
        return (0, False, 0, ((2, 'O'),)), (0, False, 0, ((1, 'H'), (1, 'H')))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Bk(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 97

    @property
    def isotopes_distribution(self):
        return FrozenDict({249: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({249: 249.074987})

    @property
    def _common_valences(self):
        return 0, 2, 3, 4

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()), (4, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Cf(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 98

    @property
    def isotopes_distribution(self):
        return FrozenDict({249: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({249: 249.074854})

    @property
    def _common_valences(self):
        return 0, 2, 3, 4

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Es(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 99

    @property
    def isotopes_distribution(self):
        return FrozenDict({252: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({252: 252.08298})

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Fm(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 100

    @property
    def isotopes_distribution(self):
        return FrozenDict({257: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({257: 257.095106})

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Md(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 101

    @property
    def isotopes_distribution(self):
        return FrozenDict({258: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({258: 258.098431})

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class No(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 102

    @property
    def isotopes_distribution(self):
        return FrozenDict({259: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({259: 259.10103})

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


class Lr(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 103

    @property
    def isotopes_distribution(self):
        return FrozenDict({266: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({266: 266.11983})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group


__all__ = ['Sc', 'Y',
           'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
           'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
