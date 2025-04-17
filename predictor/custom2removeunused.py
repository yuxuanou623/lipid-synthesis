"""SMARTS representations custom reactions."""
from synthemol.reactions.query_mol import QueryMol
from synthemol.reactions.reaction import Reaction

# To use custom chemical reactions instead of Enamine REAL reactions, replace None with a list of Reaction objects.
# If CUSTOM_REACTIONS is None, synthemol will default to the reactions in real.py.
# CUSTOM_REACTIONS: tuple[Reaction] | None = None


CUSTOM_REACTIONS = [
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[OH1][C:4]([*:5])=[O:6]')
        ],
        product=QueryMol('[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'),
        reaction_id=0
    ),
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[OH1][C:4]([*:5])=[O:6]')
        ],
        product=QueryMol('[*:5][C:4](=[O:6])[N:2]([*:1])[*:3]'),
        reaction_id=1
    ),
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[H:3]'),
            QueryMol('[*:4][N:5]([H])[*:6]')
        ],
        product=QueryMol('O=C([N:2]([*:1])[H:3])[N:5]([*:4])[*:6]'),
        reaction_id=2
    ),
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[H:3]'),
            QueryMol('[*:4][N:5]([H])[H:6]')
        ],
        product=QueryMol('O=C([N:2]([*:1])[H:3])[N:5]([*:4])[H:6]'),
        reaction_id=3
    ),
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[F,Cl,Br,I][*:4]')
        ],
        product=QueryMol('[*:1][N:2]([*:3])[*:4]'),
        reaction_id=4
    ),
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[H:3]'),
            QueryMol('[*:4][N:5]([H])[H:6]')
        ],
        product=QueryMol('O=C(C(=O)[N:2]([*:1])[H:3])[N:5]([*:4])[H:6]'),
        reaction_id=5
    ),
    Reaction(
        reactants=[
            QueryMol('[*:1][N:2]([H])[*:3]'),
            QueryMol('[*:4][N:5]([H])[H:6]')
        ],
        product=QueryMol('O=C(C(=O)[N:2]([*:1])[*:3])[N:5]([*:4])[H:6]'),
        reaction_id=7
    ),
    Reaction(
        reactants=[
            QueryMol('[OH1:1][C:2]([*:3])=[O:4]'),
            QueryMol('[F,Cl,Br,I][*:5]')
        ],
        product=QueryMol('[O:4]=[C:2]([*:3])[O:1][*:5]'),
        reaction_id=8
    ),
    # Amine + Acrylate: Michael addition of amine to α,β-unsaturated ester
    Reaction(
        reactants=[
            QueryMol('[N:1]([H])[H]'),
            QueryMol('[C:2]=[C:3][C:4](=[O:5])[O:6][*:7]')
        ],
        product=QueryMol('[N:1]([C:2]([C:3][C:4](=[O:5])[O:6][*:7]))[H]'),
        reaction_id=9
    ),
    # Amine + Epoxide: Amine opens an epoxide ring
    Reaction(
        reactants=[
            QueryMol('[N:1]([H:5])[H:6]'),  # Primary amine (-NH2) with explicit hydrogens
            QueryMol('[C:2]1[C:3][O:4]1')   # Epoxide (three-membered ring)
        ],
        product=QueryMol('[N:1]([H:5])[C:2][C:3]([O:4][H:6])'),  # β-Amino Alcohol
        reaction_id=10
    ),
    # Michael Addition (Thiol): Thiol addition to α,β-unsaturated carbonyl
    Reaction(
        reactants=[
            QueryMol('[S:1][H]'),
            QueryMol('[C:2]=[C:3][C:4](=[O:5])[O:6][*:7]')
        ],
        product=QueryMol('[S:1][C:2]([C:3][C:4](=[O:5])[O:6][*:7])'),
        reaction_id=11
    ),
    # Michael Addition (Amine-Ester + Thiol): Thiol addition to α,β-unsaturated amide
    Reaction(
        reactants=[
            QueryMol('[C:1](=[O:2])[N:3][C:4]=[C:5][C:6](=[O:7])[*:8]'),
            QueryMol('[S:9][H]')
        ],
        product=QueryMol('[C:1](=[O:2])[N:3][C:4]([S:9])[C:5]([C:6](=[O:7])[*:8])'),
        reaction_id=12
    ),
    # Amine + Carbonate: Amine attacks carbonate to form carbamate
    Reaction(
        reactants=[
        QueryMol('[N:1]([H])[H]'),
        QueryMol('[O:2][C:3](=[O:4])[O:5]')
    ],
        product=QueryMol('[N:1][C:3](=[O:4])[O:5]'),
        reaction_id=13
    ),
    # Carboxylic Acid + Alcohol → Ester (via Fischer esterification)
    Reaction(
        reactants=[
            QueryMol('[C:1](=[O:2])[OH]'),  # Carboxylic acid (-COOH)
            QueryMol('[O:3][H]')          # Alcohol (-OH with attached R group)
        ],
        product=QueryMol('[O:3][C:1](=[O:2])'),  # Correct Ester (-COO-R)
        reaction_id=15
    ),


]