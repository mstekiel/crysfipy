import numpy as np
from typing import Union

def sympy_expression_mapping():
    return {'sqrt':np.emath.sqrt, 'I':1j}


def point_group_synonyms(name: Union[str, int]):
    '''
    Return the synonyms of a given point group name.

    Contains database of Hermann-Mauguin, Schoenflies and lattice names as np.array with
    dtype = [
            ('HM', 'U8'),
            ('Sch', 'U8'),
            ('lattice', 'U16'),
            ('id', 'i4')
            ]

    Parameters:
    -----------
    name:
        Name of the space group in Hermann-Mauguin, Schoenflies or lattice type notation.
        Alternatively, integer identifyer of the point symmetry.

        
    Returns:
    --------
    Tuple of (Hermann-Mauguin, Schoenflies, lattice names, id),
    in case the name is ambiguous a list of fitting tuples.
    '''


    symmetry_synonyms = [
        ('1',   'C1',   'triclinic',   0),
        ('2',   'C2',   'monoclinic',   1),
        ('m',   'Cs',   'monoclinic',   1),
        ('2/m', 'C2h',  'monoclinic',   1),
        ('mm2', 'C2v',  'orthorhombic', 2),
        ('222', 'D2',   'orthorhombic', 2),
        ('mmm', 'D2h',  'orthorhombic', 2),
        ('4',   'C4',   'tetragonal_1', 3),
        ('-4',  'S4',   'tetragonal_1', 3),
        ('4/m', 'C4h',  'tetragonal_1', 3),
        ('422', 'D4',   'tetragonal_2', 4),
        ('4mm', 'C4v',  'tetragonal_2', 4),
        ('-42m',        'D2d',  'tetragonal_2', 4),
        ('4/mmm',       'D4h',  'tetragonal_2', 4),
        ('3',   'C3',   'trigonal_1',   5),
        ('-3',  'S6',   'trigonal_1',   5),
        ('32',  'D3',   'trigonal_2',   6),
        ('3m',  'C3v',  'trigonal_2',   6),
        ('-3m', 'D3d',  'trigonal_2',   6),
        ('6',   'C6',   'hexagonal_1',  7),
        ('-6',  'C3h',  'hexagonal_1',  7),
        ('6/m', 'C6h',  'hexagonal_1',  7),
        ('622', 'D6',   'hexagonal_2',  8),
        ('6mm', 'C6v',  'hexagonal_2',  8),
        ('-6m2',        'D3h',  'hexagonal_2',  8),
        ('6/mmm',       'D6h',  'hexagonal_2',  8),
        ('23',  'T',    'cubic_1',      9),
        ('m-3', 'Th',   'cubic_1',      9),
        ('-43m',        'Td',   'cubic_2',      10),
        ('432', 'O',    'cubic_2',      10),
        ('m-3m',        'Oh',   'cubic_2',      10)
    ]

    symmetry_db = np.array(symmetry_synonyms, dtype=[
                                        ('HM', 'U8'),
                                        ('Sch', 'U8'),
                                        ('lattice', 'U16'),
                                        ('id', 'i4')
                                        ])
    
    symmetry = None

    if isinstance(name, str):
        if name in symmetry_db['HM']:
            symmetry = symmetry_db[ symmetry_db['HM']==name ]
        if name in symmetry_db['Sch']:
            symmetry = symmetry_db[ symmetry_db['Sch']==name ]
        if name in symmetry_db['lattice']:
            symmetry = symmetry_db[ symmetry_db['lattice']==name ]
    if isinstance(name, int):
        if name in symmetry_db['id']:
            symmetry = symmetry_db[ symmetry_db['id']==name ]

    if symmetry is None:
        raise NameError(f'{name} is an unknown point group identifier.')

    return symmetry

# class CEF_tables():
# 	pass
# 	'''
# 	Hard-coded data on various CEF related objects, printing and conversion functionalities.
#
# 	'''
# 	def __init__(self):
# 		return

# 	def fill_symmetry(self):
# 		return



def composeTables(filename: str) -> str:
    '''
    Internal function that makes the tables of elements with principal values required for calculations. It is not well implemented, but it works and makes things easy.
    '''
    
    # Principal factors as implemented originally by Czech guys
    M = { 
       # ion    J     gJ
		"ce" : [2.5,  6.0/ 7.0,   -2.0/35.0              ,  2.0/315.0                  ,  0.0                            ],
        "pr" : [4.0,  4.0/ 5.0,   -2.0**2*13/3**2/5**2/11, -2.0**2/3**2/5/11**2        ,  2.0**4*17/3**4/5/7/11**2/13    ],
        "nd" : [4.5,  8.0/11.0,    7.0/1089.0            , -136.0/467181.0             , -1615.0/     42513471.0         ],
        #"pm" : [4.0,  3.0/ 5.0,    2.0*7/3/5/11**2       ,  2.0**3*7*17/3**3/5/11**3/13,  2.0**3*17*19/3**3/7/11**2/13**2], # Promethium is missing from form factors database
        "sm" : [2.5,  2.0/ 7.0,   13.0/3**2/5/7          ,  2.0*13/3**3/5/7/11         ,  0.0                            ],
        "tb" : [6.0,  3.0/ 2.0,   -1.0/99.0              ,  2.0/        16335.0        ,  1.0/(3**4*7*11**2*13)          ],
        "dy" : [7.5,  4.0/ 3.0,   -2.0/3**2/5/7          , -2.0**3/3**3/5/7/11/13      ,  2.0**2/3**3/7/11**2/13**2      ],
        "ho" : [8.0,  5.0/ 4.0,   -1.0/2/3**2/5**2       , -1.0/2/3/5/7/11/13          , -5.0/3**3/7/11**2/13**3         ],
        "er" : [7.5,  6.0/ 5.0,    4.0/(3**2*5**2*7)     ,  2.0/(3**2*5*7*11*13)       ,  8.0/(3**3*7*11**2*13)          ],
        "tm" : [6.0,  7.0/ 6.0,    1.0/3**2/11           ,  2.0**3/3**4/5/11**2	       , -5.0/3**4/7/11**2/13            ],
        "yb" : [3.5,  8.0/ 7.0,    2.0/3**2/7            , -2.0/3/5/7/11               ,  2.0**2/3**3/7/11/13            ]}
        
    # Form factors database by J.P. Brown
    jDatabase = {}
    with open(filename, 'r') as ff:
        lines = ff.readlines()
        
    for line in lines:
        if line[0] == '#':
            continue
        else:
            fields = line.replace('=',' ').split()
            name = fields[1].lower()
            keys = fields[3::2]
            values = np.array(fields[4::2], dtype=float)
            jorder = keys[0][3]
            
            if name in jDatabase:
                jDatabase[name][f'j{jorder}'] = values
                
            else:
                jDatabase[name] = {f'j{jorder}' : values}
            
    atomicDatabase = {}
    for name in M:
        if name == 'ce':
            oxidation = '2'
        else:
            oxidation = '3'
            
        extendedName = name+oxidation
        
        atomicDatabase[name] = {'J':M[name][0], 'gJ':M[name][1], 'alpha':M[name][2], 'beta':M[name][3], 'gamma':M[name][4], 'js':jDatabase[extendedName]}
     
     
    precision = 10
    out = ['atomicDatabase = {']
    for name in atomicDatabase:
        out.append(f'"{name}": {{"J":{atomicDatabase[name]["J"]}, "gJ":{atomicDatabase[name]["gJ"]}, ')
        out.append(f'\t"alpha":{atomicDatabase[name]["alpha"]:.{precision}e},\t"beta":{atomicDatabase[name]["beta"]:.{precision}e},\t"gamma":{atomicDatabase[name]["gamma"]:.{precision}e},')
        out.append(f'\t"js":{{\t"j0":{list(atomicDatabase[name]["js"]["j0"])},')
        out.append(f'\t\t"j2":{list(atomicDatabase[name]["js"]["j2"])} }}}},')
        
    out.append('}')
        
    return '\n'.join(out)
        
if __name__ == '__main__':
    pass
    # tables = composeTables('./atoms-form-factors.txt')
    # print(tables)

    # db = point_group_synonyms(6)
    # print(db)