#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
:mod:`arithmetique` module
:author: FIL - Faculté des Sciences et Technologies -  Univ. Lille <http://portail.fil.univ-lille1.fr>_
:date: décembre 2018
Fonctions arithmétiques abordées en MIAC

"""

def mon_divmod(a, b):
    '''
    :param a,b: (int) deux entiers
    :valeur renvoyée: (tuple) l'unique couple (q, r) tel que 
             a = bq + r   avec 0 <= r < bcorrespondant au quotient et au reste 
             (division euclidienne de a par b)
    :CU: a >= 0 et b > 0

    Algo : recherche du plus grand multiple de b <= a

    Exemples :

    >>> mon_divmod(10, 3)
    (3, 1)
    >>> mon_divmod(12, 3)
    (4, 0)
    '''
    q, r = 0, a
    # invariant a == b * q + r and r >= 0
    while r >= b:
        q = q + 1
        r = r - b
        # invariant a == b * q + r and r >= 0
    # a == b * q + r  and  0 <= r < a
    return (q, r)

def pgcd(a, b):
    '''
    :param a, b: (int) deux nombres entiers
    :valeur renvoyée: (int) pgcd de a et b
    :CU: (a, b) != (0, 0)
    :Exemples:

    >>> pgcd(12, 8)
    4
    >>> pgcd(0, 3)
    3
    >>> pgcd(42, 26)
    2
    >>> pgcd(-42, 26)
    2
    >>> pgcd(42, -26)
    2
    >>> pgcd(-42, -26)
    2
    >>> pgcd(47337338776803609, 114875698154581919)
    1
    '''
    a1, b1 = abs(a), abs(b)
    # invariant : pgcd(a, b) = pgcd(a1, b1)
    while b1 != 0:
        a1, b1 = b1, a1 % b1
        # invariant : pgcd(a, b) = pgcd(a1, b1)
    # pgcd(a, b) = pgcd(a1, b1) et b1 = 0
    return abs(a1)

def euclide_etendu(a, b):
    '''
    :param a, b: (int) deux nombres entiers
    :valeur renvoyée: (tuple) un triplet d'entiers (d, u, v) tels que
        - d = pgcd (a,b)
        - et d = au + bv 
    :CU: (a, b) != (0, 0) non tous deux nuls
    :Exemple:

    >>> a, b = 875667429081, 76547625781
    >>> d, u, v = euclide_etendu(a, b)
    >>> d == pgcd(a, b)
    True
    >>> d == a * u + b * v
    True
    '''
    r0, r1 = a, b
    u0, u1 = 1, 0
    v0, v1 = 0, 1
    # Invariant :   pgcd (a,b) = pgcd (r0, r1) et a x u0 + b x v0 = r0
    #               et a x u1 + b x v1 = b1
    while r1 != 0:
        q, r = divmod(r0, r1)
        r0, r1 = r1, r
        u0, u1 = u1, u0 - q*u1
        v0, v1 = v1, v0 - q*v1
        # Invariant :  pgcd (a,b) = pgcd (r0, r1) et a x u0 + b x v0 = r0 
    #  pgcd (a,b) = pgcd (a1,b1) et a x u0 + b x v0 = r0 et r1 = 0
    if r0 < 0:
        r0, u0, v0 = -r0, -u0, -v0
    return (r0, u0, v0)


def inverse_mod(a, n):
    """
    :param a, n: (int) deux entiers
    :valeur renvoyée: (int) l'inverse de a modulo n
    :CU: n non nul et a inversible modulo n 
    :Exemples:

    >>> inverse_mod(3, 7)
    5
    >>> inverse_mod(3, 26)
    9
    >>> inverse_mod(47337338776803609, 114875698154581919)
    64484911222477487
    """
    assert n != 0, 'n ne doit pas être nul'
    d, u, v = euclide_etendu(a, n)
    assert d == 1, 'a non inversible modulo n'
    return u % n


def resoud_equation(a, b, c):
    """
    :param a, b, c: (int) trois nombres entiers coefficients d'une 
               équation de la forme ax + by = c (1)
    :valeur renvoyée: (tuple) 
       - couple de deux couples (x0, y0) et (alpha, beta)
         tels que l'ensemble des solutions de l'équation (1) soit 
         l'ensemble des couples de la forme
         (x0 + alpha x k, y0 + beta x k), k entier relatif quelconque,
         si (1) a des solutions.
       - () si (1) n'a pas de solutions.
    :CU: (a, b) != (0, 0)
    :Exemples:

    >>> resoud_equation(4, 6, 5)
    ()
    >>> from random import randint
    >>> a, b = randint(1, 100), randint(-100, 100)
    >>> c = (a, b, a*b)[randint(0, 2)]
    >>> sols = resoud_equation(a, b, c)
    >>> x0, y0 = sols[0]
    >>> alpha, beta = sols[1]
    >>> all(a*(x0 + alpha*k) + b*(y0 + beta*k) - c == 0 for k in range(-20, 20))
    True
    """
    d = pgcd(a, b)
    c1, r = divmod(c, d)
    if r != 0:
        # coeff c non divisible par pgcd(a, b) => pas de solution
        return ()
    else:
        # l'équation a des solutions. On rend les coeff a et b premiers entre eux.
        a1, b1 = (a // d, b // d)
        # on cherche les coeff de Bezout (u,v) tq a1 u + b1 v = 1
        d, u, v = euclide_etendu(a1, b1)
        # on multiplie (u, v) par c1 pour obtenir une solution de a1 x0 + b1 y0 = c1
        x0, y0 = (c1 * u, c1 * v)
        return ((x0, y0), (-b1, a1))

def eratosthene(n):
    '''
    :param n: (int) 
    :return: (list)  liste des nbres premiers <= n.
    :CU: n >= 2
    :Exemples:

    >>> eratosthene(40)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    >>> len(eratosthene(100)) == 25
    True
    '''
    premiers = []
    barres = set()
    p = 2
    while p * p <= n:
        if p not in barres:
            # p est premier
            premiers.append(p)
            # barrer tous les multiples de p
            mp = p + p
            while mp <= n:
                barres.add(mp)
                mp += p
        p += 1
    for k in range (p, n + 1):
        if k not in barres:
            premiers.append(k)
    return premiers

def est_premier1(n):
    """
    :param n: (int) entier à tester
    :valeur renvoyée: (bool)
       - True si n est  un nombre premier, 
       - False sinon
    :CU: 2 <= n
    :Exemples:

    >>> est_premier1(9929)
    True
    >>> est_premier1(9927)
    False
    """
    return n in eratosthene(n)


def factorise(n):
    """
    :param n: (int) l'entier à factoriser
    :valeur renvoyée: (list) liste de couples (p, alpha) avec p nombre premier
             et alpha non nul plus grand entier tq p ^ alpha | n 
    :CU: n >= 1
    :Exemples:

    >>> factorise(1)
    []
    >>> factorise(2)
    [(2, 1)]
    >>> factorise(504)
    [(2, 3), (3, 2), (7, 1)]
    """
    m = n
    facteurs = []
    p = 2
    alpha = 0
    while m != 1:
        if m % p == 0:
            alpha += 1
            m = m // p
        else:
            if alpha != 0:
                facteurs.append((p, alpha))
            p += 1
            alpha = 0
    if alpha != 0:
        facteurs.append((p, alpha))
    return facteurs

def plus_petit_diviseur1(n):
    """
    :param n: (int) entier dont on cherche le plus petit diviseur > 1.
    :valeur renvoyée: (int) le plus petit diviseur de n (autre que 1)
    :CU: n >= 2
    :Exemples:

    >>> plus_petit_diviseur1(9927)
    3
    >>> plus_petit_diviseur1(9929)
    9929
    """
    d = 2
    while n % d != 0:
        d += 1
    return d


def plus_petit_diviseur(n):
    """
    :param n: (int) entier dont on cherche le plus petit diviseur > 1.
    :valeur renvoyée: (int) le plus petit diviseur de n (autre que 1)
    :CU: n >= 2
    :Exemples:

    >>> plus_petit_diviseur(9927)
    3
    >>> plus_petit_diviseur(9929)
    9929
    """
    d = 2
    while d * d <= n and n % d != 0:
        d += 1
    if d * d <= n:
        return d
    else:
        return n

def est_premier2(n):
    """
    :param n: (int) entier à tester
    :valeur renvoyée: (bool)
       - True si n est un nombre premier, 
       - False sinon 
    :CU: 2 <= n
    :Exemples:

    >>> est_premier2(9929)
    True
    >>> est_premier2(9927)
    False
    """
    return plus_petit_diviseur(n) == n



def expo_mod_rapide(a, b, n):
    """
    :param a: (int)
    :param b: (int) exposant
    :param n: (int) modulus
    :valeur renvoyée: (int) r dans [] tel que a^b = r (mod n)
              calcul effectué par l'exponentiation rapide.
    :CU: b >= 0, n > 0
    :Exemples:

    >>> expo_mod_rapide(2, 10, 100)
    24
    >>> expo_mod_rapide(14, 3141, 17)
    12
    """
    r, s, k = 1, a, b
    # invariant : a^b % n = r*s^k % n
    while k != 0:
        if k % 2 == 1:
            r = (r * s) % n
        s = (s * s) % n
        k = k // 2
        # invariant : a^b % n = r*s^k % n
    # a^b % n = r*s^k % n et k = 0
    return r



def est_temoin_non_primalite(a, n):
    '''
    :param a, n: (int)
    :valeur renvoyée: (bool)
      - True si a est un témoin de Fermat de non primalité de n
      - False sinon
    :CU: 1 < n 
    :Exemples:

    >>> n = 2**32 + 1
    >>> est_temoin_non_primalite(2, n)
    False
    >>> est_temoin_non_primalite(3, n)
    True
    '''
    return pgcd(a, n) == 1 and expo_mod_rapide(a, n - 1, n) != 1

def est_premier_probable_fermat(n, nbtentatives=20):
    '''
    :param n: (int) entier dont on teste la primalité (probable)
    :param nbtentatives: (int) (optionnel) nbre de candidats témoins à essayer
    :valeur renvoyée: (bool)
      - True
      - False
    :CU: n >= 3
    '''
    return not any(est_temoin_non_primalite(randrange(2, n-1), n)
                   for k in range(nbtentatives))

def est_carmichael(n):
    '''
    :param n: (int)
    :return: (bool)
       - True si n est un nombre de Carmichael
       - False sinon
    :CU: n > 0
    :Exemples:
    
    >>> any(est_carmichael(k) for k in range(1, 561))
    False
    >>> est_carmichael(561)
    True
    '''
    fact = factorise(n)
    r = len(fact)
    return (r > 1 and
            all(fact[i][1] == 1 for i in range(r)) and
            all((n - 1) % (fact[i][0] - 1) == 0 for i in range(r)))

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS, verbose=True)
    

