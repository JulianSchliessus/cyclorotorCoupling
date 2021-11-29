#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from itertools import islice
import configparser
import numpy as np
import meshio
import pathlib
from .helper import Mesh


class MBDynPrep:
    ''' A simple preprocessor that provides the data
    for the precice-mbdyn adapter. '''

    def __init__(self, case_name=None, in_mm=True):
        self.name = ''
        self.mesh = Mesh()
        self.problem_dict = dict()
        self.control_dict = dict()
        self.nodes_dict = dict()
        self.material_dict = dict()
        self.design_variable_dict = dict()
        self.main_path = str(pathlib.Path(__file__).parent.resolve())

        # TODO: case name input instead of mesh name
        if case_name:
            if os.path.dirname(case_name):
                directory = os.path.dirname(case_name)
                config_name = os.path.join(directory, 'mbdyn-config')
            else:
                exit('Wrong path')
            self.read_config(config_name)
            if 'fixed nodes' in self.nodes_dict.keys():
                self.read_gmsh(case_name, 
                               self.nodes_dict['fixed nodes'], mm_to_m=in_mm)
            else:
                self.read_gmsh(case_name, mm_to_m=in_mm)

    def read_config(self, file_name='mbdyn-config'):
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.read(file_name)
        self.control_dict = dict(dict(parser)['control'])
        self.problem_dict = dict(dict(parser)['problem'])
        self.nodes_dict = dict(dict(parser)['nodes'])
        self.material_dict = dict(dict(parser)['material'])
        self.design_variable_dict = dict(dict(parser)['Design Variables'])

    # TODO: Test changes
    def read_gmsh(self, file_name, fixed_nodes=None, mm_to_m=False):
        self.name = os.path.splitext(os.path.basename(file_name))[0]
        mesh_file = open(file_name, 'r')
        if 'MeshFormat' in mesh_file.readline():
            format_info = mesh_file.readline().split(sep=' ')
            if format_info[0] == '2.2' and format_info[1] == '0':
                print('Compatible gmsh file format found!')
                next(mesh_file)
            else:
                raise(ImportError(
                    "Mesh format of '{name}' incompatible.".format(
                        name=file_name)))
        names = list()
        nodes = None
        edges = list()
        shells = list()
        for line in mesh_file:
            if line.strip('$\n') == 'PhysicalNames':
                nlines = int(list(islice(mesh_file, 0, 1))[0])
                for block_line in islice(mesh_file, 0, nlines):
                    cells = block_line.split(sep=' ')
                    names.append((int(cells[0]), int(cells[1]),
                                  cells[2].strip('"\n')))
                # print(names)
            elif line.strip('$\n') == 'Nodes':
                nlines = int(list(islice(mesh_file, 0, 1))[0])
                nodes = np.genfromtxt(islice(mesh_file, 0, nlines),
                                      dtype=float)
                # print(nodes)
            else:
                if 'End' not in line:
                    print('Warning: Unknown Block found in mesh file.')
        mesh_file.close()

        self.mesh.name = file_name
        self.mesh.nodes = nodes[:, 1:]

        if fixed_nodes is not None:
            self.mesh.set_clamp_constraint(fixed_nodes, dead_z=False)
            
        if mm_to_m:
            self.mesh.nodes *= 0.001

    # TODO
    def read_meshio(self, file_name, fixed_nodes=None):
        meshio_mesh = meshio.read(file_name)

        assert len(meshio_mesh.points), \
            "Mesh does not contain any nodes!"
        self.mesh.nodes = meshio_mesh.points

        assert 'quad' in meshio_mesh.cells_dict.keys(), \
            "Mesh does not contain any quadrilateral elements!"
        self.mesh.shells = meshio_mesh.get_cells_type('quad')

        if 'line' in meshio_mesh.cell_dict.keys():
            self.mesh.edges = meshio_mesh.get_cells_type('line')

        if 'gmsh:physical' in meshio_mesh.cell_data_dict.keys():
            if meshio_mesh.field_data:
                edge_type = meshio_mesh.get_cell_data('gmsh:physical', 'line')
                shell_type = meshio_mesh.get_cell_data('gmsh:physical', 'quad')
                physical_names = meshio_mesh.field_data

        if fixed_nodes is not None:
            self.mesh.set_clamp_constraint(fixed_nodes)


if __name__ == '__main__':
    prep = MBDynPrep('test-files/kite-v5-1500cells.msh')
