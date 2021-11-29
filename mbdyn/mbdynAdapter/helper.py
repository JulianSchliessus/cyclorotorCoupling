#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from subprocess import Popen
import os
import glob
from mbc_py_interface import mbcNodal
from bs4 import BeautifulSoup
import precice
import numpy as np
from scipy.spatial.transform import Rotation as R
import logging
import pathlib
from ..rotation import *		#rotation.py
#import matplotlib.pyplot as plt
import math

# create logger
module_logger = logging.getLogger('adapter.helper')

class MBDynHelper:
    def __init__(self, mesh):
        self.initialized = False
        self.process = None
        self.nodal = None
        self.log_file = None
        self.main_path = str(pathlib.Path(__file__).parent.resolve())
        self.log_file_path = self.main_path + '/../../mbdyn.log'
        self.input_file_name = 'shell.mbd'
        self.mesh = mesh
        self.load_changed = False
        self.pressure = 0
        self.stresses = 0
        self.node_forces = 0
        #self.cell_forces = None
        self._debug_samples = [0]

    def _equidistant_samples(self, num_samples=6):
        self._debug_samples.append(0)
        num_nodes = self.mesh.number_of_nodes()
        interval = num_nodes/(num_samples-1)
        for i in range(1, num_samples):
            self._debug_samples.append(int(interval*i)-1)


    def initialize(self, case='shell'):
        for i in glob.glob('mbdyn/mbdynData/*.sock'):
        	os.remove(i)
        self.case = case
        self.input_file_name = self.main_path + '/../mbdynData/' + self.case + '.mbd'
        self.log_file = open(self.log_file_path, 'w')
        self.process = Popen(['mbdyn', '-f', self.input_file_name],
                             stdout=self.log_file,
                             stderr=self.log_file)
        self.process.stdin = ''
        
    def nodal_initialize(self):    
        path = self.main_path + '/../mbdynData/{name}.sock'.format(name=self.case)
        '''if len(path) > 107:
        	exit("\n\nPlease use a shorter MBDyn socket path! (Limit: 107 characters)\nCurrent path: " + path + "\nMain path:" + self.main_path + "\n\n")'''
        path = "mbdyn/mbdynData/{name}.sock".format(name=self.case)
        module_logger.debug('socket path: %s' % path)
        host = ''
        port = 0
        timeout = -1
        verbose = 1
        data_and_next = 1
        refnode = 0
        nodes = self.mesh.number_of_nodes()
        labels = 0 # 16
        rot = 256 # for rotvec 256, rotmat 512, euler 1024; see mbc.h enum MBCType
        accels = 0
        
        self.nodal = mbcNodal(path, host, port, timeout, verbose,
                              data_and_next, refnode, nodes, labels,
                              rot, accels)
    def nodal_negotiate(self):
        print("prenegotiate")
        self.nodal.negotiate()
        print("post neg")
    def nodal_negotiate_recv(self):
        self.nodal.recv()
        print("post recv")
        self.initialized = True

    def finalize(self):
        try:
            self.nodal.destroy()
            self.log_file.close()
        except AttributeError:
            print('Warning: Could not close log file or destroy mbc.')

    def get_absolute_displacement(self,absolute = False, deformation=False):
        if absolute and not deformation:
        	return self.get_nodes()
        elif deformation and not absolute:
        	print(math.sqrt((self.get_nodes()[0,0] - self.get_nodes()[1,0])**2 + (self.get_nodes()[0,1] - self.get_nodes()[1,1])**2))
        	if math.sqrt((self.get_nodes()[0,0] - self.get_nodes()[1,0])**2 + (self.get_nodes()[0,1] - self.get_nodes()[1,1])**2) > 0.02:
        		print(math.sqrt((self.get_nodes()[0,0] - self.get_nodes()[1,0])**2 + (self.get_nodes()[0,1] - self.get_nodes()[1,1])**2))
        	transform = Rotation()
        	transform.importGrid(self.get_nodes()[:11,0],self.get_nodes()[:11,1],self.get_nodes()[:11,2])
        	XA, YA, ZA = transform.rotate(angle=np.degrees(float(-self.get_rotation()[11,2])),rot_point = np.array((0,0,0)))
        	transform.importGrid(XA,YA,ZA)
        	XA, YA, ZA = transform.rotate(angle=np.degrees(float(-self.get_rotation()[5,2]+self.get_rotation()[11,2])),rot_point = np.array((0.5,0,0)))

        	mesh_displacement = np.array((XA.flatten(),YA.flatten(),ZA.flatten())).T - self.mesh.nodes[:11,:]
        	
        	'''self.transform.importGrid(self.mbdyn.get_nodes()[i*npp:npp-1+i*npp,0],self.mbdyn.get_nodes()[i*npp:npp-1+i*npp,1],self.mbdyn.get_nodes()[i*npp:npp-1+i*npp,2])
XA, YA, ZA = self.transform.rotate(angle=np.degrees(float(-self.mbdyn.get_rotation()[npp-1+i*npp,2])),rot_point = np.array((0,0,0)))
self.transform.importGrid(XA,YA,ZA)
XA, YA, ZA = self.transform.rotate(angle=np.degrees(float(-self.mbdyn.get_rotation()[5+i*npp,2]+self.mbdyn.get_rotation()[npp-1+i*npp,2])),rot_point = self.mbdyn.mesh.nodes[npp-1+i*npp,:])
self.displacement.append((np.array((XA.flatten(),YA.flatten(),ZA.flatten())).T.copy() - self.mbdyn.mesh.nodes[i*npp:npp-1+i*npp,:]).copy())'''
        	
        	return mesh_displacement
        elif absolute and deformation:
        	exit("wrong condition")
        else:
        	return self.get_nodes() - self.mesh.nodes

    def get_nodes(self):
        if self.initialized:
            return np.reshape(self.nodal.n_x, (-1, 3))
        else:
            return self.mesh.nodes

    def get_forces(self):
        return np.reshape(self.nodal.n_f, (-1, 3))
        
    def get_rotation(self):
    	return np.reshape(self.nodal.n_theta, (-1, 3))

    def set_forces(self, forces):
        self.node_forces = forces
        #import ipdb; ipdb.set_trace()
        self.nodal.n_f[:] = np.ravel(forces)

    #TODO: create option for stresses
    def set_pressure(self, pressure):
        self.pressure = float(pressure)
        self.load_changed = True

    #TODO: something about the greyed out part breaks things
    def calc_pressure_forces(self, forces=0, relaxation=1, limiting=10):

        module_logger.debug('rotvec from mbdyn: \n %s' % self.nodal.n_theta)

        node_normals_weighted = np.zeros((self.mesh.number_of_nodes(), 3))

        pressure_forces = node_normals_weighted * self.pressure

        if not isinstance(limiting, type(None)):
            max_value_pressure = np.max(np.linalg.norm(pressure_forces, axis=1))
            if max_value_pressure > limiting:
                pressure_forces = self.node_forces
            if not isinstance(self.node_forces, (int, float)):
                max_value_fluid = np.max(np.linalg.norm(forces, axis=1))
                if max_value_fluid > limiting:
                    forces = forces / max_value_fluid * 0.3

        if relaxation != 1:
            new_forces = self.node_forces + \
                (pressure_forces + forces - self.node_forces) * relaxation
        else:
            new_forces = pressure_forces + forces

        # self.node_forces = np.multiply(fix_force, self.pressure)
        forces_norm = np.linalg.norm(new_forces, axis=1)
        module_logger.debug(
            'min, max, sum forces after pressure applied:\n{}, {}, {}'.format(
                np.min(forces_norm), np.max(forces_norm),
                np.sum(forces_norm)))
        module_logger.debug(
            'forces after pressure applied sample:\n{}'.format(
                new_forces[self._debug_samples,:]))

        self.set_forces(new_forces)


    def solve(self, converged=False):
        if self.nodal.send(converged):
            module_logger.debug('on send')
            return True
        if self.nodal.recv():
            module_logger.debug('on recv')
            return True
        return False

    #TODO
    def solve_static(self, tolerance=1e-6, max_iterations=10000,
                     write=True):
        previous_position = 0
        for i in range(max_iterations):
            self.calc_pressure_forces(relaxation=0.3, limiting=20000)
            if self.solve(True):
                return True
            current_position = self.get_absolute_displacement()
            two_norm_diff = np.linalg.norm(
                current_position - previous_position)
            previous_position = current_position
            module_logger.debug('Finished iteration: {}/{}, displacement two-norm diff: {}/{}'.format(
                i, max_iterations, two_norm_diff, tolerance))
            if two_norm_diff < tolerance and i > 500:
                print('Converged in {}/{} iterations'.format(
                    i, max_iterations))
                return True
        print('No convergence in {} iterations'.format(max_iterations))
        return False

    def solve_initial(self, tolerance=5e-6, max_iterations=10000,
                     write=True):

        previous_position = 0

        # calculate static force magnitude on each node
        self.calc_pressure_forces()
        node_forces_mag = np.linalg.norm(self.node_forces, axis=1)[:, np.newaxis]

        # calculate node normals
        def node_normals(xyz=self.get_nodes()):
            normals = np.zeros((self.mesh.number_of_nodes(), 3))
            return normalize_vectors(normals)

        # calculate static force in new direction
        def new_force():
            return node_normals() * node_forces_mag

        for i in range(max_iterations):
            if self.solve(True):
                return True

            current_position = self.get_absolute_displacement()

            two_norm_diff = np.linalg.norm(
                current_position - previous_position)

            previous_position = current_position

            module_logger.debug('Finished iteration: {}/{}, displacement two-norm diff: {}/{}'.format(
                i, max_iterations, two_norm_diff, tolerance))

            update = new_force()
            update *= ((i+1)/200) if i < 200 else 1
            self.node_forces = update
            self.set_forces(update)

            forces_norm = np.linalg.norm(update, axis=1)
            module_logger.debug(
                'min, max, sum forces after pressure applied:\n{}, {}, {}'.format(
                    np.min(forces_norm), np.max(forces_norm),
                    np.sum(forces_norm)))
            module_logger.debug(
                'forces after pressure applied sample:\n{}'.format(
                    update[self._debug_samples,:]))

            if two_norm_diff < tolerance and i > 500:
                module_logger.debug('Converged in {}/{} iterations'.format(
                    i, max_iterations))
                return True

        module_logger.debug('No convergence in {} iterations'.format(max_iterations))

        return False

    # TODO
    def get_node_normals(self):
        normals = np.zeros((self.mesh.number_of_nodes(),3))
        normals = normalize_vectors(normals)

        return normals


class PreciceHelper:
    def __init__(self, path):
        self.interface = None
        self.config_path = path
        self.dimensions = 0
        self.num_vertices = 0
        self.vertex_ids = 0
        self.quad_ids = 0
        self.displacement_id = 0
        self.displacement = None
        self.rotation_id = 0
        self.rotation = None
        self.force_id = 0
        self.force = None
        self.time_step = 0

    def setup_interface(self, solver_name='Solid'):
        print(solver_name, self.config_path)
        self.solver_name = solver_name
        self.interface = precice.Interface(
            solver_name, str(self.config_path), 0, 1)

    def configure_interface(self, nodes, grid_name='Solid-Mesh',
                            quads=None):
        self.num_vertices = len(nodes)
        self.dimensions = self.interface.get_dimensions()

        mesh_id = self.interface.get_mesh_id(grid_name)
        vertices = nodes

        self.displacement = np.zeros((self.num_vertices, self.dimensions))
        self.force = np.zeros((self.num_vertices, self.dimensions))

        self.vertex_ids = self.interface.set_mesh_vertices(
            mesh_id, vertices)
        module_logger.debug('precice vertex ids:\n %s' % str(self.vertex_ids))

        if not isinstance(quads, type(None)):
            for ids in quads:
                self.quad_ids = self.interface.set_mesh_quad_with_edges(
                    mesh_id, ids[0], ids[1], ids[2], ids[3])

        self.displacement_id = self.interface.get_data_id(
            'Displacement_Data', mesh_id)
        self.rotation_id = self.interface.get_data_id(
            'Rotation_Data', mesh_id)
        self.force_id = self.interface.get_data_id(
            'Force_Data', mesh_id)

        self.time_step = self.interface.initialize()

        if self.interface.is_read_data_available():
        	if self.solver_name== 'Structure_Solver':
        		self.interface.read_block_vector_data(self.force_id,
        			self.vertex_ids)
        	else:
        		exit("Error")
        		self.interface.read_block_vector_data(self.displacement_id,
        			self.vertex_ids)

    def initialize_data(self):
        if self.interface.is_action_required(
                precice.action_write_initial_data()):
            if self.solver_name== 'Structure_Solver':
            	self.interface.write_block_vector_data(self.displacement_id,
            		self.vertex_ids,
            		self.displacement)
            	self.interface.write_block_vector_data(self.rotation_id,
            		self.vertex_ids,
            		self.rotation)
            else:
            	exit("Error")
            	self.interface.write_block_vector_data(self.force_id,
            		self.vertex_ids,
            		self.force)
            self.interface.mark_action_fulfilled(
                precice.action_write_initial_data())

        self.interface.initialize_data()

    def get_participant_name_from_xml(self):
        with open(self.config_path, 'r', encoding='utf8') as file:
            content = file.read()
            soup = BeautifulSoup(content, 'xml')
            participants = soup.find_all("participant")
            for solver in participants:
                if 'Force' in solver.find('read-data').attrs['name']:
                    name = solver.attrs['name']
        return name

    def get_mesh_name_from_xml(self):
        with open(self.config_path, 'r', encoding='utf8') as file:
            content = file.read()
            soup = BeautifulSoup(content, 'xml')
            participants = soup.find_all("participant")
            for solver in participants:
                if 'Force' in solver.find('read-data').attrs['name']:
                    mesh_names = solver.find_all('use-mesh')
                    for mesh in mesh_names:
                        if 'provide' in mesh.attrs:
                            name = mesh.attrs['name']
        return name

    def advance_time(self):
        print("MBDyn Adapter: Advancing in time")
        self.time_step = self.interface.advance(self.time_step)

    def read_data(self):
        if self.interface.is_read_data_available():
        	print("Reading DATA")
        	if self.solver_name== 'Structure_Solver':
        		self.force = self.interface.read_block_vector_data(
        			self.force_id, self.vertex_ids)
        	else:
        		exit("Error")
        		self.displacement = self.interface.read_block_vector_data(
            		self.displacement_id, self.vertex_ids)
        else:
         	print("No readable DATA")

    def write_data(self, write_data_displacement,write_data_rotation):
        if self.interface.is_write_data_required(self.time_step):
        	if self.solver_name== 'Structure_Solver':
        		self.interface.write_block_vector_data(
        			self.displacement_id, self.vertex_ids, write_data_displacement)
        		self.interface.write_block_vector_data(
        			self.rotation_id, self.vertex_ids, write_data_rotation)
        	else:
        		exit("Error")
        		self.interface.write_block_vector_data(
        			self.force_id, self.vertex_ids, write_data)

class Mesh:
    def __init__(self):
        self.name = ''
        self.nodes = np.array(None)
        self.node_constraints = np.array(None)
        self.node_orientations = np.array(None)
        self.edges = np.array(None)
        self.edge_names = np.array(None)
        self.shells = np.array(None)
        self.shell_names = np.array(None)


    def constraints_from_edge_names(self):
        self.node_constraints = np.full((len(self.nodes), 6), False,
                                        dtype='?')
        for idx in range(len(self.edges)):
            cur_constraints = np.full(6, False, dtype='?')
            cur_name = self.edge_names[idx].casefold()
            if 'fix' in cur_name:
                cur_name = cur_name.replace('fix', '')
                if cur_name == 'all':
                    cur_constraints[:] = True
                else:
                    if 'x' in cur_name:
                        cur_constraints[0] = True
                    if 'y' in cur_name:
                        cur_constraints[1] = True
                    if 'z' in cur_name:
                        cur_constraints[2] = True
                    if 'a' in cur_name:
                        cur_constraints[3] = True
                    if 'b' in cur_name:
                        cur_constraints[4] = True
                    if 'c' in cur_name:
                        cur_constraints[5] = True
            for node in self.edges[idx]:
                self.node_constraints[node, :] += cur_constraints

    def number_of_nodes(self):
        return len(self.nodes)

    def calc_node_orientation(self, unchanged='z', clean_unchanged=True,
                              return_normal=False, flip_normal=1):
        valid_unchanged = ['x', 'y', 'z']
        if unchanged not in valid_unchanged:
            raise ValueError('unchanged must be one of {}.'.format(
                unchanged))


        if clean_unchanged:
            cell_normals[:, valid_unchanged.index(unchanged)] = 0

        cell_normals = normalize_vectors(cell_normals)

        node_normals = np.zeros((self.number_of_nodes(), 3))

        node_normals = normalize_vectors(node_normals)

        # TODO: which coordinate points outwards
        global_frame = np.zeros((2, 3))
        global_frame[0, 1] = 1
        global_frame[1, 2] = 1

        local_frame = global_frame.copy()

        orientation = np.zeros((self.number_of_nodes(), 3))

        for i, normal in enumerate(node_normals):
            local_frame[0, :] = normal
            rotation = R.align_vectors(global_frame, local_frame)
            orientation[i, :] = rotation[0].as_euler('xyz')

        self.node_orientations = orientation

        if return_normal:
            return orientation, node_normals
        return orientation

    #TODO: fix if, test it
    def set_clamp_constraint(self, fixed_nodes, dead_z=False):
        assert(isinstance(fixed_nodes, (slice, list, int)))
        if not self.node_constraints.any():
            self.node_constraints = np.full(
                (self.number_of_nodes(), 6),  False, dtype='?')
        self.node_constraints[fixed_nodes, :3] = True
        if dead_z:
            self.node_constraints[:, 2] = True

def normalize_vectors(vectors):
    length = np.linalg.norm(vectors, axis=1)
    return np.divide(vectors.transpose(), length).transpose()


if __name__ == "__main__":
    pass
