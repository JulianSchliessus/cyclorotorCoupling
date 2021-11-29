import numpy as np
#import matplotlib.pyplot as plt

class Rotation:
	def makeGrid(self,Lx,Ly,Lz,meshDensity):
		self.Lx_min, self.Lx_max = Lx
		self.Ly_min, self.Ly_max = Ly
		self.Lz_min, self.Lz_max = Ly
		
		exit("deprecated!")
		### Not tested!!! -> change to z coordinate necessary
		self.YA, self.XA = np.mgrid[self.Ly_min:self.Ly_max:(meshDensity*1j), self.Lx_min:self.Lx_max:(meshDensity*1j)]
		### Not tested!!!
		
		return np.array((self.XA,self.YA,self.ZA))
	
	def importGrid(self,XA,YA,ZA):
		self.XA = XA
		self.YA = YA
		self.ZA = ZA
		self.Lx_min = np.amin(self.XA)
		self.Lx_max = np.amax(self.XA)
		self.Ly_min = np.amin(self.YA)
		self.Ly_max = np.amax(self.YA)
		self.Lz_min = np.amin(self.ZA)
		self.Lz_max = np.amax(self.ZA)
		#return np.array((self.Lx_min,self.Lx_max,self.Ly_min,self.Ly_max,self.Lz_min,self.Lz_max))
	
	def centerGrid(self,x0,y0,z0):
		
		#translate
		self.XA = self.XA - x0
		self.YA = self.YA - y0
		self.ZA = self.ZA - z0

		return np.array((self.XA,self.YA,self.ZA))
		
	def translateBack(self,x0,y0,z0):

		#translate
		self.XAprim = self.XAprim + x0
		self.YAprim = self.YAprim + y0
		self.ZAprim = self.ZAprim + z0

		return np.array((self.XAprim,self.YAprim,self.ZAprim))
	
	def rotate(self, rot_point = None, angle = None):
		#about z-axis
		
		
		#center of gravity if no center entered
		x0 = rot_point[0] if (rot_point is not None) else (self.Lx_max+self.Lx_min)/2.
		y0 = rot_point[1] if (rot_point is not None) else (self.Ly_max+self.Ly_min)/2.
		z0 = rot_point[2] if (rot_point is not None) else (self.Lz_max+self.Lz_min)/2.
		#import ipdb; ipdb.set_trace()
		self.centerGrid(x0,y0,z0)
		
		angle = angle if (angle is not None) else 45.
		ang_rot = np.radians(angle)

		#about z-axis
		#counterclockwise
		RotMatrix = np.array([[np.cos(ang_rot),  -np.sin(ang_rot),0],
                      		[np.sin(ang_rot), np.cos(ang_rot),0],
                      		[0,0,1]])
		#print(str(np.degrees(ang_rot))+"°")
		#rotate

		self.XAprim, self.YAprim, self.ZAprim = np.einsum('ji, mni -> jmn', RotMatrix, np.dstack([self.XA, self.YA, self.ZA]))
		#import ipdb; ipdb.set_trace()
		self.XA, self.YA, self.ZA = self.translateBack(x0,y0,z0)
		return np.array((self.XA, self.YA,self.ZA))
