import copy
import random

from base_complex.singularity_structure import SingularityStructure

# build a component mesh based on the hybrid mesh
# can a component mesh equal to the hybrid mesh?


# new version
class HybridSingularityGraph(SingularityStructure):
    def __init__(self, mesh, iteration_num=0):
        SingularityStructure.__init__(self, mesh, iteration_num)

    # hybrid singularity graph should be about to support component elements
