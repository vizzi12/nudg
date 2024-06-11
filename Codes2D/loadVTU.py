# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 09:40:01 2024

@author: victo
"""

import meshio
import numpy as np
import scipy


for i in range(4):
    mesh = meshio.read("SWE_000"+str(i)+"_0000.vtu")
    mesh.point_data['tri']=mesh.cells_dict['triangle']
    
    mesh.point_data['points']=mesh.points[:,0:2]
    
    scipy.io.savemat('test'+str(i)+'.mat',mesh.point_data)
    print(str(i))



"""
for i in range(2):    
    mesh = meshio.read("SWECircularDam1_000"+str(i)+"_0001.vtu")
    mesh.point_data['tri']=mesh.cells_dict['triangle']
    
    mesh.point_data['points']=mesh.points[:,0:2]
    
    scipy.io.savemat('circularTest'+str(i)+'.mat',mesh.point_data)
    print(str(i))
    """