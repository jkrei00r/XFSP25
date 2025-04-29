# region imports
import math
import random as rnd
from copy import deepcopy as dc


# endregion

# region class definitions
class Position:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Position(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Position(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar):
        return Position(self.x * scalar, self.y * scalar, self.z * scalar)

    def __truediv__(self, scalar):
        return Position(self.x / scalar, self.y / scalar, self.z / scalar)

    def mag(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        mag = self.mag()
        if mag > 0:
            self.x /= mag
            self.y /= mag
            self.z /= mag

    def random_direction(self):
        vec = Position(rnd.gauss(0, 1), rnd.gauss(0, 1), rnd.gauss(0, 1))
        vec.normalize()
        return vec


class Molecule:
    def __init__(self, position=Position(), weight=14):
        self.position = position
        self.weight = weight


class Macromolecule:
    def __init__(self, target_N=1000):
        self.N = max(1, int(rnd.gauss(target_N, 0.1 * target_N)))
        self.segment_length = 25.4e-9  # Adjusted to match example scale
        self.mers = []
        self.com = Position()
        self.e2e = 0.0
        self.rog = 0.0
        self.weight = 14 * self.N

    def build_chain(self):
        current_pos = Position()
        self.mers = [Molecule(position=dc(current_pos))]

        for _ in range(self.N - 1):
            direction = Position().random_direction()
            current_pos += direction * self.segment_length
            self.mers.append(Molecule(position=dc(current_pos)))

        # Calculate metrics
        self.calculate_com()
        self.e2e = (self.mers[0].position - self.mers[-1].position).mag()
        self.calculate_rog()

    def calculate_com(self):
        total_weight = sum(m.weight for m in self.mers)
        com = Position()
        for m in self.mers:
            com += m.position * m.weight
        self.com = com / total_weight

    def calculate_rog(self):
        total_weight = sum(m.weight for m in self.mers)
        sum_sq = sum(m.weight * (m.position - self.com).mag() ** 2 for m in self.mers)
        self.rog = math.sqrt(sum_sq / total_weight)
# endregion
