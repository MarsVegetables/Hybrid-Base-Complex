import copy
import random

from base_complex.singularity_structure import SingularityStructure


# new version
class HybridSingularityGraph(SingularityStructure):
    def __init__(self, mesh, iteration_num=0):
        SingularityStructure.__init__(self, mesh, iteration_num)

    # hybrid singularity graph should be able to support component elements
