from odbAccess import *
import numpy as np

odb = openOdb(path='j1.odb')

# here is the data at day 8 
lastFrame = odb.steps['migration1'].frames[-1]


CellDensity1 = lastFrame.fieldOutputs['NT11']
CellDensity2 = lastFrame.fieldOutputs['NT12']
CellDensity3 = lastFrame.fieldOutputs['NT13']
Disp = lastFrame.fieldOutputs['U']



sideline = odb.rootAssembly.instances['PART-1-1'].nodeSets['SIDELINE']
n1 = odb.rootAssembly.instances['PART-1-1'].nodeSets['N1']

cell1= CellDensity1.getSubset(region=sideline)
cell2= CellDensity2.getSubset(region=sideline)
cell3= CellDensity3.getSubset(region=sideline)
disp_n1= Disp.getSubset(region=n1)


dispFile = open('NT11_1.txt','w')
for c in cell1.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('NT12_1.txt','w')
for c in cell2.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('NT13_1.txt','w')
for c in cell3.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('disp_1.txt','w')
for c in disp_n1.values:
    dispFile.write('%10.4E\n' % (c.data[0]))
dispFile.close() 



# here is the data at day 17  
lastFrame = odb.steps['migration2'].frames[-1]


CellDensity1 = lastFrame.fieldOutputs['NT11']
CellDensity2 = lastFrame.fieldOutputs['NT12']
CellDensity3 = lastFrame.fieldOutputs['NT13']
Disp = lastFrame.fieldOutputs['U']



sideline = odb.rootAssembly.instances['PART-1-1'].nodeSets['SIDELINE']
n1 = odb.rootAssembly.instances['PART-1-1'].nodeSets['N1']

cell1= CellDensity1.getSubset(region=sideline)
cell2= CellDensity2.getSubset(region=sideline)
cell3= CellDensity3.getSubset(region=sideline)
disp_n1= Disp.getSubset(region=n1)



dispFile = open('NT11_2.txt','w')
for c in cell1.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('NT12_2.txt','w')
for c in cell2.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('NT13_2.txt','w')
for c in cell3.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('disp_2.txt','w')
for c in disp_n1.values:
    dispFile.write('%10.4E\n' % (c.data[0]))
dispFile.close() 


# here is the data at day 17  
lastFrame = odb.steps['migration3'].frames[-1]


CellDensity1 = lastFrame.fieldOutputs['NT11']
CellDensity2 = lastFrame.fieldOutputs['NT12']
CellDensity3 = lastFrame.fieldOutputs['NT13']
Disp = lastFrame.fieldOutputs['U']



sideline = odb.rootAssembly.instances['PART-1-1'].nodeSets['SIDELINE']
n1 = odb.rootAssembly.instances['PART-1-1'].nodeSets['N1']

cell1= CellDensity1.getSubset(region=sideline)
cell2= CellDensity2.getSubset(region=sideline)
cell3= CellDensity3.getSubset(region=sideline)
disp_n1= Disp.getSubset(region=n1)



dispFile = open('NT11_3.txt','w')
for c in cell1.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('NT12_3.txt','w')
for c in cell2.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('NT13_3.txt','w')
for c in cell3.values:
    dispFile.write('%10.4E\n' % (c.data))
dispFile.close()

dispFile = open('disp_3.txt','w')
for c in disp_n1.values:
    dispFile.write('%10.4E\n' % (c.data[0]))
dispFile.close() 







