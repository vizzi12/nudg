# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:40:45 2024

@author: victo
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 09:40:01 2024

@author: victo
"""

import meshio
import numpy as np
import scipy



for i in range(1,37):
    mesh = meshio.read("SWE"+str(i)+"_0000_0001.vtu")
    mesh.point_data['tri']=mesh.cells_dict['triangle']
    
    mesh.point_data['points']=mesh.points[:,0:2]
    
    scipy.io.savemat('test'+str(i)+'.mat',mesh.point_data)
    print(str(i))



