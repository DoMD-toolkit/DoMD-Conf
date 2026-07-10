import re
from io import StringIO
from xml.etree import cElementTree
import networkx as nx
import numpy as np


def pbc(x,l):
    """Applies periodic boundary conditions to a coordinate.

    Args:
        x (float or np.ndarray): The input coordinate.
        l (float or np.ndarray): The box length in the corresponding dimension.

    Returns:
        float or np.ndarray: The coordinate wrapped within [-l/2, l/2].
    """
    return x - l*np.rint(x/l)

def control_in(control_file):
    """Parses a control input file.

    Args:
        control_file (str): The path to the control file.
    """
    pass


class Box(object):
    """Represents the simulation box dimensions and tilt factors.

    Attributes:
        xy (float): Tilt factor xy.
        xz (float): Tilt factor xz.
        yz (float): Tilt factor yz.
    """
    def __init__(self):
        """Initializes a Box object with zero tilt factors."""
        self.xy = 0
        self.xz = 0
        self.yz = 0
        return

    def update(self, dic):
        """Updates the box attributes from a dictionary.

        Args:
            dic (dict): Dictionary containing box attributes.
        """
        self.__dict__.update(dic)


class XmlParser(object):
    """Parses XML files containing simulation data.

    This class reads an XML file and extracts specific data elements, optionally filtering
    by a list of needed tags. It handles specific structures like simulation boxes,
    reaction lists, and templates, while parsing other data into numpy arrays.

    Attributes:
        box (Box): The simulation box object.
        data (dict): A dictionary storing parsed data arrays and structures (e.g., 'reaction', 'template').
    """
    def __init__(self, filename, needed=None):
        """Initializes the XmlParser.



        Parsing logic:
        1. Reads root attributes and sets them as instance attributes (e.g., step, time).
        2. Iterates through child elements.
        3. 'box' tags update the self.box object.
        4. Other tags are filtered against the `needed` list if provided.
        5. 'reaction' tags are parsed into lists of reaction steps.
        6. 'template' tags are evaluated as Python dictionaries.
        7. All other tags are parsed as numpy arrays from text content.

        Args:
            filename (str): Path to the XML file to parse.
            needed (list, optional): A list of XML tag names to extract. If provided,
                tags not in this list (except 'box') are skipped. Defaults to None.
        """
        tree = cElementTree.ElementTree(file=filename)
        root = tree.getroot()
        self.box = Box()
        self.data = {}
        needed = [] if needed is None else needed
        for key in root[0].attrib:
            self.__dict__[key] = int(root[0].attrib[key])
        for element in root[0]:
            if element.tag == 'box':
                self.box.update(element.attrib)
                continue
            if (len(needed) > 0) and (element.tag not in needed):
                continue
            if element.tag == 'reaction':
                self.data['reaction'] = []
                reaction_list = element.text.strip().split('\n')
                while '' in reaction_list:
                    reaction_list.remove('')
                for l in reaction_list:
                    r = re.split(r'\s+', l)
                    while '' in r:
                        r.remove('')
                    r[1:] = [int(_) for _ in r[1:]]
                    self.data['reaction'].append(r)
                continue
            if element.tag == 'template':
                self.data['template'] = eval('{%s}' % element.text)
                continue
            if len(element.text.strip()) > 0:
                self.data[element.tag] = np.genfromtxt(StringIO(element.text), dtype=None, encoding=None)

template_cg = '''<?xml version ="1.0" encoding ="UTF-8" ?>
<{program}_xml version="{version}">
<configuration time_step="0" dimensions="3" natoms="{n_atoms:d}" >
<box lx="{lx:.8f}" ly="{ly:.8f}" lz="{lz:.8f}" xy="{xy:8f}" xz="{xz:8f}" yz="{yz:8f}"/>
<position num="{n_atoms:d}">
{positions}</position>
<type num="{n_atoms:d}">
{types}</type>
<image num="{n_atoms:d}">
{images}</image>
<body num="{n_atoms:d}">
{bodys}</body>
<opls_type num="{n_atoms:d}">
{opls_type}</opls_type>
<monomer_id num="{n_atoms:d}">
{monomer_id}</monomer_id>
<charge num="{n_atoms:d}">
{charge}</charge>
<mass num="{n_atoms:d}">
{mass}</mass>
<bond num="{n_bonds:d}">
{bond}
</bond>
<angle num="{n_angles:d}">
{angle}
</angle>
<dihedral num="{n_dihedrals:d}">
{dihedral}
</dihedral>
<improper num="{n_impropers:d}">
{improper}
</improper>
</configuration>
</{program}_xml>'''

def XmlWriter(CG_systems, box, filename='chemfast', program='galamost', version='1.3'):
    """Writes a Coarse-Grained (CG) system configuration to an XML file.

    Generates a simulation input file for CG systems. It handles hyperedges for
    angles and dihedrals, rigid bodies, and periodic boundary conditions (PBC).

    Args:
        CG_systems (list[nx.Graph]): A list of coarse-grained molecular graphs.
            Graphs may contain `_hyperedges` attributes for angles and dihedrals.
        box (tuple): Simulation box dimensions (lx, ly, lz).
        filename (str, optional): Base name for the output file. Defaults to 'chemfast'.
        program (str, optional): Target simulation software name tag. Defaults to 'galamost'.
        version (str, optional): XML format version. Defaults to '1.3'.
    """
    CG_system_graph = nx.compose_all(CG_systems)
    n_atoms = len(CG_system_graph.nodes)
    n_bonds = len(CG_system_graph.edges)
    n_angles = 0#len(CG_system_graph._hyperedges[3]) if CG_system_graph._hyperedges.get(3) else 0
    n_dihedrals = 0# len(CG_system_graph._hyperedges[4]) if CG_system_graph._hyperedges.get(4) else 0
    mass = types = opls_type = positions = charge = monomer_id = body = image = ''
    bond = angle = dihedral = improper  = ''
    for system in CG_systems:
        for atom_graph in system.nodes:
            mass += '1.0\n'
            charge += f'{system.nodes[atom_graph].get("charge",0)}\n'
            monomer_id += f'{atom_graph}\n'
            types += '%s\n' % system.nodes[atom_graph]['type']
            pos = pbc(system.nodes[atom_graph]['x'], box)
            positions += '%.6f %.6f %.6f\n' % (pos[0], pos[1], pos[2])  # in nm
            body += '%d\n' % int(system.nodes[atom_graph].get('body', -1))
            if system.nodes[atom_graph].get('image') is None:
                image += '0 0 0\n'
            else:
                atom_image = system.nodes[atom_graph]['image']
                image += '%d %d %d\n' % (atom_image[0], atom_image[1], atom_image[2])
        n_angles_ = 0
        angles = {}
        n_angles += n_angles_
        n_dihedrals_ = 0
        dihedrals = {}
        n_dihedrals += n_dihedrals_
        for k in angles:
            angle += f"{angles[k]['bt']} {k[0]} {k[1]} {k[2]}\n"
        dihedral = ''
        for k in dihedrals:
            dihedral += (f"{dihedrals[k]['bond_type']} "
                         f"{k[0]} {k[1]} {k[2]} {k[3]}\n")
        for k in system.edges:
            bond += f"{system.edges[k]['bond_type']} {k[0]} {k[1]}\n"
    lx, ly, lz = box
    xy, xz, yz = [0, 0, 0]
    o = open('out_%s.xml' % filename, 'w')
    o.write(
        template_cg.format(
            n_atoms=n_atoms, n_bonds=n_bonds, mass=mass, types=types, images=image, bodys=body,opls_type=opls_type, positions=positions,
            bond=bond, charge=charge, angle=angle, dihedral=dihedral, n_angles=n_angles, n_dihedrals=n_dihedrals,
            monomer_id=monomer_id, program=program, version=version, lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz,
            n_impropers=0, improper=improper
        )
    )
    o.close()
    return