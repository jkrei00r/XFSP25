# polymerClasses.py

"""
This module simulates a freely jointed polymer chain (e.g., poly(ethylene)) in 3D space.

Key Features:
-------------
- Implements random-walk-based generation of polymer chains using a fixed bond length.
- Computes physical metrics including:
    - Center of mass (COM)
    - End-to-end distance (Re)
    - Radius of gyration (Rg)
- Designed for educational or lightweight simulation purposes.

Notes on Units:
---------------
In previous versions, the segment (bond) length was set to 200 nm, which is several orders
of magnitude too long for polyethylene. The corrected value used here is 1.54 nm, matching
the typical C–C bond length. With this change, output metrics (Rg and Re) will scale
realistically for polymers with N ≈ 1000 monomers.

"""

import math
import random as rnd
from copy import deepcopy as dc


class Position:
    """
    A 3D vector to represent positions or directions in space.
    Supports vector arithmetic and normalization.
    """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        """Vector addition"""
        return Position(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        """Vector subtraction"""
        return Position(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar):
        """Scalar multiplication"""
        return Position(self.x * scalar, self.y * scalar, self.z * scalar)

    def __truediv__(self, scalar):
        """Scalar division"""
        return Position(self.x / scalar, self.y / scalar, self.z / scalar)

    def mag(self):
        """Returns the Euclidean magnitude of the vector"""
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        """
        Converts this vector to a unit vector (in-place).
        Does nothing if the magnitude is zero.
        """
        mag = self.mag()
        if mag > 0:
            self.x /= mag
            self.y /= mag
            self.z /= mag

    def random_direction(self):
        """
        Returns a new unit vector pointing in a random 3D direction.
        Uses a Gaussian distribution for isotropy.
        """
        vec = Position(rnd.gauss(0, 1), rnd.gauss(0, 1), rnd.gauss(0, 1))
        vec.normalize()
        return vec


class Molecule:
    """
    Represents a monomer in the polymer chain with position and weight.

    Attributes:
    -----------
    position : Position
        The 3D coordinates of the monomer.
    weight : float
        Mass or relative weight of the monomer (default = 14 g/mol for CH2 group).
    """

    def __init__(self, position=Position(), weight=14):
        self.position = position
        self.weight = weight


class Macromolecule:
    """
    Represents a full polymer chain made up of monomers.

    Attributes:
    -----------
    N : int
        Number of monomer units in the chain.
    segment_length : float
        Fixed bond length between monomers (1.54 nm = 1.54e-9 m).
    mers : list of Molecule
        The ordered list of monomers (molecules) forming the chain.
    com : Position
        Center of mass of the chain.
    e2e : float
        End-to-end distance in meters.
    rog : float
        Radius of gyration in meters.
    weight : float
        Total molecular weight of the polymer.
    """

    def __init__(self, target_N=1000):
        # Allow minor variation in chain length via normal distribution
        self.N = max(1, int(round(rnd.gauss(target_N, 0.1 * target_N))))

        # Use realistic C–C bond length for poly(ethylene): 1.54 nm
        self.segment_length = 1.54e-9

        self.mers = []
        self.com = Position()
        self.e2e = 0.0
        self.rog = 0.0
        self.weight = 14 * self.N  # Approx. CH2 repeat unit

    def build_chain(self):
        """
        Constructs the polymer chain as a random walk in 3D space.
        Each monomer is separated by a fixed bond length in a random direction.
        """
        current_pos = Position()
        self.mers = [Molecule(position=dc(current_pos))]

        # Add each subsequent monomer based on a random direction
        for _ in range(self.N - 1):
            direction = Position().random_direction()
            current_pos += direction * self.segment_length
            self.mers.append(Molecule(position=dc(current_pos)))

        # Calculate structural metrics
        self.calculate_com()
        self.e2e = (self.mers[-1].position - self.mers[0].position).mag()
        self.calculate_rog()

    def calculate_com(self):
        """
        Computes the center of mass of the polymer chain, weighted by monomer mass.
        """
        total_weight = sum(m.weight for m in self.mers)
        com = Position()
        for m in self.mers:
            com += m.position * m.weight
        self.com = com / total_weight

    def calculate_rog(self):
        """
        Computes the radius of gyration of the polymer chain.
        Defined as the RMS distance of monomers from the center of mass.
        """
        total_weight = sum(m.weight for m in self.mers)
        sum_sq = sum(m.weight * (m.position - self.com).mag() ** 2 for m in self.mers)
        self.rog = math.sqrt(sum_sq / total_weight)


