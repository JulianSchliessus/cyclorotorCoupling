import logging
#from matplotlib import pyplot as plt
import precice
import time
import subprocess
from .rotation import *		#rotation.py
from .csvreader import *	#csvreader.py
import pathlib
import numpy as np
import decimal
from .mbdynAdapter.prep import MBDynPrep
import socket
import errno
from scipy import interpolate
import colorama
import struct
from .socketTools import *# SocketTools as s

class MbdynAdapter:

	#initialize values
	colorama.init(autoreset=True)
	
	from .initialize import __init__
	
	#run preCICE coupling
	def run(self):
		s = SocketTools()
		interface = precice.Interface(self.participant_name, self.configuration_file_name,
				                      self.solver_process_index, self.solver_process_size)

		for i in self.mesh_name:
			self.mesh_id.append(interface.get_mesh_id(str(i)))
		dimensions = interface.get_dimensions()
		
		for i in range(self.patches):
			self.patch.append(Rotation())
			XA, YA, ZA = np.hsplit(csvImport(self.current_path + '/../' + self.fluid_folder_name + '/patch' + str(i+1) + '.csv'),3)	#if mesh is in file
			self.patch[i].importGrid(XA,YA,ZA)
			self.vertices.append(np.array([XA.flatten(),YA.flatten(),ZA.flatten()]).T)
			
		if not (len(self.vertices) == len(self.mesh_id)):
			exit("Cannot match mesh to vertices")
		
		#interpolate from beam to airfoil
		for i in range(self.patches):
			x = self.vertices[i][:,0]
			x = np.roll(x, -1)
			x = np.split(x,2)
			y = self.vertices[i][:,1]
			y = np.roll(y, -1)
			y = np.split(y,2)
			y_min = np.array((np.array(y[0]).min(axis=0),np.array(y[1]).min(axis=0))).max(axis=0)
			y_max = np.array((np.array(y[0]).max(axis=0),np.array(y[1]).max(axis=0))).min(axis=0)
			y_center = y_min + (y_max-y_min)/2
			num_of_points = 11
			
			#x = f(y)
			f = [interpolate.interp1d(y[0],x[0]),interpolate.interp1d(y[1],x[1])]
			ynew = [np.linspace(y_min,y_max,11),np.linspace(y_min,y_max,11)]
			xnew = [f[0](ynew[0]),f[1](ynew[1])]
			
			#x[0] rechts
			#x[1] links
			#von links nach rechts
			airfoil_mesh = np.concatenate(np.array((np.array((xnew[1],ynew[1],np.zeros((xnew[1].size)))) ,np.array((xnew[0],ynew[0],np.zeros((xnew[0].size)))))),axis=1)
			self.vertices[i] = airfoil_mesh.T

			self.vertex_ids.append(interface.set_mesh_vertices(self.mesh_id[i],self.vertices[i]))
			self.read_data_id.append(interface.get_data_id(self.read_data_name[i], self.mesh_id[i]))
			self.write_data_id.append(interface.get_data_id(self.write_data_name[i], self.mesh_id[i]))
		
		self.mbdyn.initialize(case=self.mbdyn_prep.name)
		
		print("trying to connect to socket stream...")
		#MBDyn drive caller socket
		self.tvc_socket = s.driveConnect()
		self.exchange_socket = s.patchConnect(self.patches+1)
		print(u'\u2705' + " " + colorama.Style.BRIGHT + colorama.Back.WHITE + colorama.Fore.GREEN + str(self.patches+2) + " socket(s) connected")
		
		dt = interface.initialize()
		
		self.input.update_time_step(dt)
		
		print("preInit")
		self.mbdyn.nodal_initialize()
		#send tvc socket data
		s_out_bufsize = 8
		i = 0
		#tvc =  self.iteration * 0.02
		tvc = 0
		buf_out = bytearray(struct.pack("d",tvc))
		self.tvc_socket[i].sendall(buf_out, s_out_bufsize)
		print("tvc: " + str(tvc))
		#send tvc socket data
		i = 0
		#tvc =  self.iteration * 0.02
		tvc = 0
		buf_out = bytearray(struct.pack("d",tvc))
		self.tvc_socket[i].sendall(buf_out, s_out_bufsize)
		print("tvc: " + str(tvc))
		self.mbdyn.nodal_negotiate()
		#TODO: Is first iteration step here?
		self.mbdyn.nodal_negotiate_recv()


		self.transform = Rotation()
		
		self.current_time_step = decimal.Decimal(str(dt))

		while interface.is_coupling_ongoing():
			
			if interface.is_action_required(
				    precice.action_write_iteration_checkpoint()):
				print("DUMMY: Writing iteration checkpoint")
				interface.mark_action_fulfilled(
				    precice.action_write_iteration_checkpoint())

			self.read_data.clear()
			self.force_tensor.clear()
			for i in range(self.patches):
				if interface.is_read_data_available():
					self.read_data.append(interface.read_block_vector_data(self.read_data_id[i], self.vertex_ids[i]))
				elif float(self.current_time_step) == float(decimal.Decimal(str(dt))):
					self.read_data.append(self.vertices[i]*0)
				else:
					exit("no data provided!")
				self.force_tensor.append(np.split(self.read_data[i],2)[0]*0) # prepare force_tensor for entries
				self.force_tensor[i] = np.concatenate((self.force_tensor[i],self.force_tensor[i][:2,:]*0)) # prepare force_tensor for entries, add empty anchor & ground force node
				for k,j in enumerate(np.split(self.read_data[i],2)):
					self.force_tensor[i] += np.concatenate((j,j[:2,:]*0)) # fill values into force_tensor, add empty anchor & ground force node
					
			#print(np.vstack(self.force_tensor))
			#import ipdb; ipdb.set_trace()
			self.mbdyn.set_forces(np.vstack(self.force_tensor))
			
			target_tvc_angle = 0
			tvc_speed = 1*10**(-1.8)
			if float(self.current_time_step) >= 0.25 and float(self.current_time_step) < 0.5:
				target_tvc_angle = np.pi/4
			elif float(self.current_time_step) >= 0.5 and float(self.current_time_step) < 0.75:
				target_tvc_angle = 0
			elif float(self.current_time_step) >= 0.75 and float(self.current_time_step) < 1:
				target_tvc_angle = -np.pi/4
			elif float(self.current_time_step) >= 1 and float(self.current_time_step) < 1.25:
				target_tvc_angle = np.pi/2
			elif float(self.current_time_step) >= 1.25:
				target_tvc_angle = -np.pi/2
			
			if target_tvc_angle > tvc and not target_tvc_angle < tvc + tvc_speed:
				tvc = tvc + tvc_speed
			elif target_tvc_angle < tvc and not target_tvc_angle > tvc - tvc_speed:
				tvc = tvc - tvc_speed
			#import ipdb; ipdb.set_trace()
			#tvc = self.iteration * 0.02
			
			#send tvc socket data
			i = 0
			buf_out = bytearray(struct.pack("d",tvc))
			self.tvc_socket[i].sendall(buf_out, s_out_bufsize)
			print("\ntvc: " + str(tvc))
			print("\ntarget tvc angle: " + str(target_tvc_angle))
			
			if self.mbdyn.solve(False):
				self.module_logger.debug('Something went wrong!')
				#MBDyn diverges
				break
				
			#displacement_real = self.mbdyn.get_absolute_displacement()
			#displacement = self.mbdyn.get_absolute_displacement(False,True)
			#import ipdb; ipdb.set_trace()
			#nodes per patch
			npp = self.mbdyn.mesh.number_of_nodes()/self.patches
			if npp % 2 != 0 and npp % 1 != 0:
				exit("please check number of nodes / patches!")
			npp = int(npp)
			#import ipdb; ipdb.set_trace()
			self.displacement_absolute.clear()
			self.displacement.clear()
			self.rotation.clear()
			self.message.clear()
			self.centre_of_gravity.clear()
			
			i=0
			self.displacement_absolute.append(self.mbdyn.get_absolute_displacement()[i*npp:(i+1)*npp,:].copy())
			self.rotation.append(self.mbdyn.get_rotation()[i*npp:(i+1)*npp,:].copy())
			
			self.message.append("{:.4f}".format(self.displacement_absolute[i][npp-1+i*npp,0]) + ";")
			self.message[i] += "{:.4f}".format(self.displacement_absolute[i][npp-1+i*npp,1]) + ";"
			self.message[i] += "{:.4f}".format(self.displacement_absolute[i][npp-1+i*npp,2]) + ";"
			self.message[i] += "{:.4f}".format(self.rotation[i][npp-1+i*npp,0]*0) + ";"
			self.message[i] += "{:.4f}".format(self.rotation[i][npp-1+i*npp,1]*0) + ";"

			#self.message[i] += "{:.4f}".format(tvc)
			
			#no Farfield rotation
			self.message[i] += "{:.4f}".format(tvc*0)
			
			self.displacement_absolute.clear()
			self.rotation.clear()
			for i in range(self.patches):
				
				#displacement
				#TODO: change reference system of rotation -> to ground_2 / ground_3
				#
				self.displacement_absolute.append(self.mbdyn.get_absolute_displacement()[i*npp:(i+1)*npp,:].copy())
				
				#TODO: change rot Point? No -> its static
				#transform rotor rotation of most deformed element
				#TODO: move CofG of mbyn.get_nodes, delete changes at rot center below
				#import ipdb; ipdb.set_trace()
				self.transform.importGrid(
				self.mbdyn.get_nodes()[i*npp:npp-2+i*npp,0]-self.mbdyn.get_nodes()[npp-1+i*npp,0]+self.mbdyn.mesh.nodes[npp-1+i*npp,0],
				self.mbdyn.get_nodes()[i*npp:npp-2+i*npp,1]-self.mbdyn.get_nodes()[npp-1+i*npp,1]+self.mbdyn.mesh.nodes[npp-1+i*npp,1],
				self.mbdyn.get_nodes()[i*npp:npp-2+i*npp,2]-self.mbdyn.get_nodes()[npp-1+i*npp,2]+self.mbdyn.mesh.nodes[npp-1+i*npp,2],)
				XA, YA, ZA = self.transform.rotate(angle=np.degrees(float(-self.mbdyn.get_rotation()[npp-2+i*npp,2])),rot_point = np.array((0,0,0)))
				self.transform.importGrid(XA,YA,ZA)
				#TODO: why npp-2 in rotation? ... already corrected for CG in message
				#-> delta between deformation and oscillation
				XA, YA, ZA = self.transform.rotate(angle=np.degrees(float(-self.mbdyn.get_rotation()[5+i*npp,2]+self.mbdyn.get_rotation()[npp-2+i*npp,2])),rot_point = self.mbdyn.mesh.nodes[npp-2+i*npp,:])
				self.displacement.append((np.array((XA.flatten(),YA.flatten(),ZA.flatten())).T.copy() - self.mbdyn.mesh.nodes[i*npp:npp-2+i*npp,:]).copy())

				
				# centre of gravity
				#import ipdb; ipdb.set_trace()
				# i == 0 (rechts)
				# i == 1 (links)
				#import ipdb; ipdb.set_trace()
				if i == 0:
					self.centre_of_gravity.append(np.array((0.5,0,0)))
				else:
					self.centre_of_gravity.append(np.array((-0.5,0,0)))
				#self.centre_of_gravity.append(self.mbdyn.get_nodes()[(i+1)*npp-1])
				#import ipdb; ipdb.set_trace()
				### Rotation um z achse: -pi ... pi
				self.rotation.append(self.mbdyn.get_rotation()[i*npp:(i+1)*npp,:].copy())
				#import ipdb; ipdb.set_trace()
				
				#displacement is independent of CofG
				#rotation must be just oscillation #TODO
				#TODO: Change CofG to Anchor point 12
				#TODO: Adjust MBDyn to Paraview centre
				
				self.message.append("{:.4f}".format(self.displacement_absolute[i][5,0]) + ";")
				self.message[i+1] += "{:.4f}".format(self.displacement_absolute[i][5,1]) + ";"
				self.message[i+1] += "{:.4f}".format(self.displacement_absolute[i][5,2]) + ";"
				self.message[i+1] += "{:.4f}".format(self.rotation[i][5,0]*0) + ";"
				self.message[i+1] += "{:.4f}".format(self.rotation[i][5,1]*0) + ";"
				self.message[i+1] += "{:.4f}".format(self.mbdyn.get_rotation()[5+i*npp,2]-self.mbdyn.get_rotation()[npp-1+i*npp,2])
				
				'''self.message.append(str(float(0)) + ";")
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0)) + ";"
				self.message[i] += str(float(0))'''
			
			for i in range(self.patches+1):
				print('Sending')
				print(self.message[i])
				self.exchange_socket[i].sendall(str.encode(self.message[i]))
				print('Sent')
				data = self.exchange_socket[i].recv(150)
				if str(data) == "b''":
					self.exchange_socket[i].close()
					exit("Client exit")
				split_string = str(data).split("\\n")
				data = split_string[0][2:]
				if str(data) == "exit":
					for j in range(self.patches):
						self.exchange_socket[j].close()
					exit("Client exit")
				print('Received', repr(data))
	
			self.write_data.clear()
			for i in range(self.patches):
				#mesh_displacement = np.zeros((num_of_points,2))
				mesh_displacement = np.concatenate(np.array((self.displacement[i][:11,:],self.displacement[i][:11,:])),axis=0)
				#TODO: preCICE issue: undefined motion
				#self.write_data.append(mesh_displacement)
				self.write_data_tmp = np.array((mesh_displacement))
				#import ipdb;ipdb.set_trace()
				#transform coordinate system (because of the farfield rotation)
				self.transform.importGrid(self.write_data_tmp[:,0],self.write_data_tmp[:,1],self.write_data_tmp[:,2])
				XA, YA, ZA = self.transform.rotate(angle=np.degrees(float(self.rotation[i][5,2])))
				print("rotation angle (oscillation):")
				print(np.degrees(float(self.rotation[i][5,2])))
				print("rotation angle (rotor):")
				print(np.degrees(float(self.rotation[i][11,2])))
				
				self.write_data.append((np.array((XA.flatten(),YA.flatten(),ZA.flatten())).T).copy())
				
				
				#No preCICE
				#self.write_data.append((np.array((XA.flatten()*0,YA.flatten()*0,ZA.flatten()*0)).T).copy())

			
			if interface.is_write_data_required(dt):
				for i in range(self.patches):
					interface.write_block_vector_data(self.write_data_id[i], self.vertex_ids[i], self.write_data[i])
			
			print("DUMMY: Advancing in time")
			dt = interface.advance(dt)
			self.current_time_step += decimal.Decimal(str(dt))
			if interface.is_action_required(
				    precice.action_read_iteration_checkpoint()):
				print("DUMMY: Reading iteration checkpoint")
				interface.mark_action_fulfilled(
				    precice.action_read_iteration_checkpoint())
				exit("ERROR MBDyn got no values")
			else:
				#previous displacement = displacement.copy()
				
				#MBDyn advance
				if self.mbdyn.solve(True):
					print("diverged")
					break

			self.iteration += 1
		print("preCICE finalizing")
		interface.finalize()
		print("DUMMY: Closing python solver dummy...")
