from odbAccess import *
import csv

odb = openOdb('Job-5.odb')

step = odb.steps['Step-1']
frame = step.frames[-1]

# 节点集
instance = odb.rootAssembly.instances['PART-1-1']
nodeset = instance.nodeSets['BOTTOM']

# RF3
rf = frame.fieldOutputs['RF'].getSubset(region=nodeset)

# 建立RF字典
rfDict = {}
for value in rf.values:
    rfDict[value.nodeLabel] = value.data[2]   # RF3

# 输出CSV
with open('RF3_coordinate.csv','w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Node','X','Y','RF3'])

    for node in nodeset.nodes:
        x,y,z = node.coordinates
        rf3 = rfDict.get(node.label,0.0)
        writer.writerow([node.label,x,y,rf3])

odb.close()