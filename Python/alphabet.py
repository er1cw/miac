#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
:mod:`alphabet` module
:author: FIL - Faculté des Sciences et Technologies -  Univ. Lille <http://portail.fil.univ-lille1.fr>_
:date: mai 2018
:dernière modification: mars 2019

définition d'une classe Alphabet

"""

class AlphabetError(Exception):
    """
    Exception used by methods

    * ``__init__``
    * ``__getitem__``
    * ̀`index``

    of class :class:`Alphabet`.
    """
    def __init__(self, msg):
        self.message = msg


class Alphabet(object):
    """
    Creates an alphabet object.
    - if alpha is None, alphabet has 26 latin capital letters A, B, ..., Z
    - else alpha must be a string of (at least two) distinct characters

    >>> latin = Alphabet()
    >>> len(latin)
    26
    >>> ''.join(latin[i] for i in range(len(latin)))
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    >>> latin.index('Z')
    25
    >>> decimal = Alphabet('0123456789')
    >>> len(decimal)
    10
    >>> decimal[3]
    '3'
    >>> set(decimal) == set('0123456789')
    True
    >>> greek = Alphabet(''.join(chr(0x03B1 + i) for i in range(25)))
    >>> print(greek)
    αβγδεζηθικλμνξοπρςστυφχψω
    """

    def __init__(self, alpha=None):
        if not alpha is None:
            if isinstance(alpha, str):
                n = len(alpha)
                if n < 2:
                    raise AlphabetError('alphabet must be of length at least 2')
                if n != len(set(alpha)):
                    raise AlphabetError('characters in alphabet must be distinct')
                self.__alphabet__ = alpha
                self.__size__ = n
            else:
                raise AlphabetError('alphabet must be a string')
        else:
            self.__alphabet__ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            self.__size__ = 26
            
    def __repr__(self):
        return self.__alphabet__

    def __len__(self):
        return self.__size__

    def __getitem__(self, i):
        """
        :param i: an index
        :type i: int
        :return: character of index i in self
        :UC: 0 <= i < len(self)
        :raise: :class:`AlphabetError` if i is not an integer or out of index range
        """
        try:
            return self.__alphabet__[i]
        except IndexError:
            raise AlphabetError('index out of range')
        except TypeError:
            raise AlphabetError('alphabet indices must be integers')

    def __iter__(self):
        alpha = Alphabet(str(self))
        alpha.__index__ = 0
        return alpha

    def __next__(self):
        if self.__index__ == len(self):
            raise StopIteration
        char = self[self.__index__]
        self.__index__ += 1
        return char
    
    def index(self, char):
        """
        :param char: 
        :type: str
        :return:
        :rtype: int
        :CU: char must be a string of length 1, raise AlphabetError otherwise
        """
        try:
            assert isinstance(char, str) and len(char) == 1
            return self.__alphabet__.index(char)
        except AssertionError:
            raise AlphabetError('letter to find must be a string of length 1')
        except ValueError:
            raise AlphabetError('letter {:s} not in alphabet'.format(char))

## Dictionnaire d'alphabets prédéfinis
ALPHABETS = {
    'CAPITAL_LATIN' : Alphabet(),
    'CAPITAL_LATIN_SPACE' : Alphabet(str(Alphabet()) + ' '),
    'LOWER_LATIN' : Alphabet(str(Alphabet()).lower()),
    'DECIMAL_DIGITS' : Alphabet('0123456789'),
    'BINARY_DIGITS' : Alphabet('01'),
    'HEXADECIMAL_DIGITS' : Alphabet('0123456789ABCDEF'),
    'ALPHANUM' : Alphabet(str(Alphabet()) + '0123456789'),
    'PRINTABLE_ASCII' : Alphabet(''.join([chr(k) for k in range(32, 127)])),
    'ASCII' : Alphabet(''.join([chr(k) for k in range(128)])),
    'LOWER_GREEK' : Alphabet(''.join(chr(0x03B1 + i) for i in range(25)))
}


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS, verbose=False)
    
    
