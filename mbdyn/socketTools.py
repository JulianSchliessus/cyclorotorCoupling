import time
import socket

class SocketTools:
	def driveConnect(self):
		#MBDyn drive caller
		exchange_socket = []
		number_of_tvc = 1
		connection_process = False
		attempt_failed = False
		i = 0
		while not (connection_process and i == number_of_tvc):
			exchange_socket.append(socket.socket(socket.AF_UNIX, socket.SOCK_STREAM))
			try:
				exchange_socket[i].connect("mbdyn/mbdynData/drive{:02d}.sock".format(i))
				connection_process = True
				i += 1
			except:
				pass
			time.sleep(0.1)
		return exchange_socket
		
	def patchConnect(self,patches):
		#OpenFOAM socket connection
		exchange_socket = []
		connection_process = False
		attempt_failed = False
		i = 0
		# i == 0 is background mesh
		while not (connection_process and i == patches):
			exchange_socket.append(socket.socket(socket.AF_UNIX, socket.SOCK_STREAM))
			try:
				exchange_socket[i].connect(("./exchange" + str(i) + ".sock"))
				connection_process = True
				i += 1
			except:
				pass
			time.sleep(0.1)
		return exchange_socket
