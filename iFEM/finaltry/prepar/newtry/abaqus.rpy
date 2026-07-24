# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2025 replay file
# Internal Version: 2024_09_20-21.00.46 RELr427 198590
# Run by 13104 on Fri Jul 24 11:27:10 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.16602, 1.16667), width=171.637, 
    height=115.733)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('./export_nodal_u_ur.py', __main__.__dict__)
#: Model: D:/test/iFEM/finaltry/prepar/newtry/job-4.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          1
#: Number of Steps:              1
#: ========== ODB ???????? ==========
#: ODB:D:\test\iFEM\finaltry\prepar\newtry\job-4.odb
#: ????Step-1
#: ????-1
#: ???1.000000000000e+00
#: ???PART-1-1
#: ???????157
#: U ??????157
#: U ???6??????0
#: ?? UR ??????157
#: ?????157
#: ?????D:\test\iFEM\finaltry\prepar\newtry\Abaqus_All_Nodal_U_UR.csv
#: ==========================================
print('RT script done')
#: RT script done
