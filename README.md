# map_touchup
Add back difference density to improve bulk solvent.

Inspired by the PLATON squeeze package (REF), map_touchup adds back a probability-weighted version of the mFo-DFc map into the bulk solvent. The probability weight is derived by passing the sigma height of each voxel through the error function raised to the power of the number Shannon voxels in the map. This estimates the probability that at least one feature of such height occurred by chance. Prominent features are given a weigh of 1, low-lying peaks 0, and peaks near protein atoms are also suppressed. Initially, a high B factor is applied to the mFo-DFc map, followed by refinement with the updated solvent. This cycle is then repeated with a smaller B factor. This approach suppresses noise and prevents low-angle features from pushing noise peaks above the inclusion threshold. In all but a few cases so far application of map_touchup improves Rfree and geometry scores, which is remarkable because free-flagged hkls were not included in the mFo-DFc map.
