# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from functools import wraps


class cached_property:
    """
    A property that is only computed once per instance and then replaces itself
    with an ordinary attribute. Deleting the attribute resets the property.
    Source: https://github.com/bottlepy/bottle/commit/fa7733e075da0d790d809aa3d2f53071897e6f76
    """

    def __init__(self, func):
        self.__doc__ = getattr(func, "__doc__")
        self.func = func

    def __get__(self, obj, cls):
        if obj is None:
            return self

        value = obj.__dict__[self.func.__name__] = self.func(obj)
        return value


def cached_method(func):
    name = f'_cached_method_{func.__name__}'

    @wraps(func)
    def wrapper(self):
        try:
            return self.__dict__[name]
        except KeyError:
            value = self.__dict__[name] = func(self)
            return value
    return wrapper


def cached_args_method(func):
    name = f'_cached_args_method_{func.__name__}'

    @wraps(func)
    def wrapper(self, *args):
        try:
            cache = self.__dict__[name]
        except KeyError:
            value = func(self, *args)
            self.__dict__[name] = {args: value}
            return value
        try:
            return cache[args]
        except KeyError:
            value = cache[args] = func(self, *args)
            return value
    return wrapper


__all__ = ['cached_property', 'cached_method', 'cached_args_method']
