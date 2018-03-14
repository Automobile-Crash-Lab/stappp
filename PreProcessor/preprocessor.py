#coding:utf-8
#Abaqus information
#The elements of the bridge: T3D2, S4R, C3D8R, B31
#----------输入文件位置------------
# node 加上 node[2] - 0 正常节点 - 1 为fixed boundary的点 其中自由度要全部固定[1,1,1,0,0,0]
# - 2 为 tie的点 这个时候需要控制其全局坐标 node[3]
import sys
modecontrol = 1  # 控制半带宽和稀疏矩阵求解的参数
# modecontrol = 1 正常求解 2 打开半带宽求解 3 打开稀疏矩阵求解 
path = sys.path[0]
filename = 'Merge-44'      # Here we need to know the fliename of the abaqus inp file
location = path[2:].replace('\\','/') + '/' + filename + '.inp'
read = open(location,'r')

#-----------读入文件到line列表-----------
#node: 在python中，文件最好全部一次性读入内存这样效率较高
#读取文件
readline = 1    #读入变量（成功读入则继续while语句）
line = []       #文件内容

while readline:
    readline=read.readline()
    line.append(readline[:-1])  #删除最后的回车符号
read.close
line_len = len(line)            #INP文件全长度

#-----------读取特征关键词位置-------
#关键词位置信息
keylocation = [0 for x in range(line_len)] 
#注意，keylocation中为所有行的关键词种类的定义信息，最后的关键词位置信息会存储在keyloc中
pos_endallinstance = 0
for i in range(0,line_len-1):
    if line[i].find('** PART INSTANCE')!=-1:           
        keylocation[i] = 1
    if line[i].find('*Node')!=-1:           
        keylocation[i] = 2
    if line[i].find('*Element')!=-1:        
        keylocation[i] = 3
    if line[i].find('*Material')!=-1:
        keylocation[i] = 5
    if line[i].find('*Boundary')!=-1:
        keylocation[i] = 6
    if line[i].find('*Instance')!=-1:
        keylocation[i] = 7
    if line[i].find('*End Instance')!=-1:
        keylocation[i] = 8
        pos_endallinstance = i
    #if line[i].find('*Nset')!=-1:       #注意只从 最后一个instance后定义的nset开始计数 需要调整
     #   keylocation[i] = 9              #需要仿照Element的方式计算出nest的结束位置
    if line[i].find('*Tie')!=-1:        #*Tie只有两行：标题和Tie的nset所以只需要记录其开始位置
        keylocation[i] = 11

for i in range(pos_endallinstance,line_len-1):       # pos_endallinstance给出instance结束的位置也给出了
    if line[i].find('*Nset')!=-1:                   # 我们需要开始记录的nset的位置
        keylocation[i] = 9 

j = 0
for i in range(0,line_len-1):
    if keylocation[i] == 3:
        j = 1
        continue
    if j == 1:
        if line[i+1][:1] == '*':
            keylocation[i] = 4            #每个element结束位置--4
            j = 0

j = 0
for i in range(0,line_len-1):
    if keylocation[i] == 9:
        j = 1
        continue
    if j == 1:
        if line[i+1][:1] == '*':
            keylocation[i] = 10            #每个nset结束位置--10
            j = 0
       
keynum = max(keylocation)    
keyloc = [[] for x in range(keynum)]          
for i in range(0,line_len-1):
    for j in range(0,keynum):
        if keylocation[i] == (j+1):
            keyloc[j].append(i)
            
partnum = len(keyloc[0])     #PART个数

#print('keyloc',keyloc)

# keyloc[0]:Part  keyloc[1]:Node keyloc[2]:Element keyloc[3]:Element结束位置
# keyloc[4]:Material keyloc[5]:Boundary  keyloc[6]:Instance  keyloc[7]:Instance结束
# keyloc[8]:Nset   keyloc[9]:Tie        

#-----------读取固定边界信息-------------
def readboundary(b_set):
    i = 0
    j = 0
    boundary = []   
    for i in range(line_len):
        if line[i].find('*Nset, nset='+b_set)!=-1:
            j = 1
            continue
        if j==1:
            temp = line[i]
            while temp.find(',')!=-1:
                boundary.append(int(temp[:temp.find(',')]))
                temp = temp[temp.find(',')+1:]
            boundary.append(int(temp))
            if line[i+1].find('*')!=-1:
                j = 0
    return boundary

# keyloc[5]中存储的就是bounary的信息读取出边界的nset 再去读取节点号
boundary = []
b_line = line[keyloc[5][0]+1]
b_set = b_line[:b_line.find(',')]
boundary=readboundary(b_set)
#print('boundary',boundary)    
        
#-----------读取NODE信息-------------
def readnode(beg,end):
    
    #节点信息
    node_num = []   #NODE信息-node序号
    node_xyz = []   #NODE信息-节点xyz坐标
    #临时变量
    i = 0
    j = 0
    
    for i in range(beg,end+1):
    
        nodeline = line[i][line[i].find(',')+1:]        #'          -5.,          -5.,          10.' in 8H
        node_num.append(int(line[i][:line[i].find(',')]))
        node_xyz.append([])                             #为了能添加内容所初始化
        
        while nodeline.find(',') != -1 :
            loc = nodeline.find(',')                    #找出','的位置
            node_xyz[j].append(float(nodeline[:loc]))
            nodeline = nodeline[loc+1:]
        node_xyz[j].append(float(nodeline[:]))
        j = j + 1
    
    return node_num,node_xyz
    
    
#-----------读取ELEMENT信息-------------
def readelement(beg,end):
    
    #单元信息
    element_num = []    #ELEMENT信息-element序号
    element_name = []   #ELEMENT信息-element名称
    element_node = []   #ELEMENT信息-单元节点号
    #临时变量
    i = 0
    j = 0
    
    for i in range(beg,end+1):
    
        elementline = line[i][line[i].find(',')+1:]    #' 10, 11, 14, 13,  1,  2,  5,  4' in 8H
        element_num.append(int(line[i][:line[i].find(',')]))
        element_node.append([])                        #为了能添加内容所初始化
        
        while elementline.find(',') != -1 :
            loc = elementline.find(',')                #找出','的位置
            element_node[j].append(int(elementline[:loc]))
            elementline = elementline[loc+1:]
        element_node[j].append(int(elementline[:]))
        j = j + 1
    
    return element_num,element_node
    

#-----------读取PART名称-------------
def readpartname(beg):
    
    partname = line[beg][line[beg].find(':'):]
    return partname


#-----------读取ELEMENT名称----------
def readelementname(beg):
    
    elementname = line[beg][line[beg].find('type')+5:]
    return elementname
    
    
#-----------读取MATERIAL信息----------
def readmaterial(beg):
    
    material = [[] for x in range(3)] 
    density = line[beg+2]
    elastic = line[beg+4]
    material[0] = line[beg][line[beg].find('name')+5:]
    material[1] = float(density[:density.find(',')])
    material[2].append(float(elastic[:elastic.find(',')]))
    material[2].append(float(elastic[elastic.find(',')+1:]))
    return material
    
#-----------读取PART信息-------------

#部件信息     
partname = readpartname(keyloc[0][0])

node = readnode(keyloc[1][0]+1,keyloc[2][0]-1)
#print('node',node[1][0])
#node[0][i] 第i个节点的全局编号
#node[1][i] 第i个节点的全局坐标

node_b = []         #Node_Boundary
for i in range(len(node[0])):
    if node[0][i] in boundary:
        node_b.append([1,1,1,1,1,1])   #固定全是1
    else:
        node_b.append([0,0,0,0,0,0])   #自由的话先给定前面的自由度自由为0
        #[0,1,0,1,0,1][1,1,0,0,0,1]
        
element = [[] for x in range(len(keyloc[2]))]
for i in range(0,len(keyloc[2])):
    element[i].append(readelementname(keyloc[2][i]))
    element[i].append(readelement(keyloc[2][i]+1,keyloc[3][i]))

#print('element',element[2])
#element[i][0] i= 0,1,2,3 如 i=1 'B31' 单元的种类
#element[i][1] 单元的信息
    #element[i][1][0] element的序号
    #element[i][1][1] element的坐标

material = []
for i in range(0,len(keyloc[4])):
    material.append(readmaterial(keyloc[4][i]))

#print('material',material[1][2])
#material[i][0] 材料的名字
#material[i][1] 材料的密度
#material[i][2] 材料的力学参数list 杨氏模量和泊松比 

#-----------找到8H单元两种属性的分别适用的点-------------
#在INP建模中 对于C3D8R单元其被Pier和Riverbank同时使用，但是其在材料的Material属性
#上不同所以需要找出使用不同材料属性的分界点
num_append1=[]
num_append2=[]
for i in range(len(node[1])):
    if node[1][i] == [0,10,0]:
        node1 = i+1     #节点1
        num_append1.append(i)
    elif node[1][i] == [0,-10,0]:
        node2 = i+1     #节点2
        num_append2.append(i)

node1_1 = len(node[0])+1
node2_2 = node1_1 + 1       

node[0].append(node1_1)
node[1].append(node[1][node1-1])
node_b.append([0,0,0,1,1,1])
node[0].append(node2_2)
node[1].append(node[1][node2-1])
node_b.append([0,0,0,1,1,1])
    

for i in range(len(keyloc[2])):
    if element[i][0] == 'C3D8R':
        for j in range(len(element[i][1][0])):
            for k in range(len(element[i][1][1][j])):
                if element[i][1][1][j][k]==node1:
                    element[i][1][1][j][k] = node1_1
                elif element[i][1][1][j][k]==node2:
                    element[i][1][1][j][k] = node2_2    
                     
#------读取按照PART区分的Merge-i文件得到两种8H单元的单元分界信息--------
filename = 'Merge-4'
location = path[2:].replace('\\','/') + '/' + filename + '.inp'
read = open(location,'r')

#读取文件
readline = 1    #读入变量（成功读入则继续while语句）
line_i = []       #文件内

#读入文件到line
while readline:
    readline=read.readline()
    line_i.append(readline[:-1])  #删除最后的回车符号
    i = i+1
read.close
line_len = len(line_i)            #INP文件全长度
#print(len(line))

C3D8_pier = 0
for i in range(line_len):
    if line_i[i].find('** PART INSTANCE: PART-PIER-1')!=-1: 
        j = 1
        continue
    if j==1 and line_i[i].find('*Nset')!=-1:
        C3D8_pier = int(line_i[i-1][:line_i[i-1].find(',')])
        break

#-----------输出材料对应单元信息-------------
material_e = [['floor'],['Pier'],['SupportBeam'],['RiverBank'],['Cables']] 
for i in range(len(keyloc[4])):
    if material[i][0] == 'STEEL':
        material_e[4].append(material[i])
        for i in range(len(keyloc[2])):
            if element[i][0] == 'T3D2':
                material_e[4].append([element[i][1][0][0],element[i][1][0][-1]])
    elif material[i][0] == 'ALUMINUM':
        material_e[2].append(material[i])
        for i in range(len(keyloc[2])):
            if element[i][0] == 'B31':
                material_e[2].append([element[i][1][0][0],element[i][1][0][-1]]) 
    elif material[i][0] == 'GRANITE':
        material_e[3].append(material[i])
        for i in range(len(keyloc[2])):
            if element[i][0] == 'C3D8R':
                material_e[3].append([C3D8_pier+1,element[i][1][0][-1]]) 
    elif material[i][0] == 'CONCRETE':
        material_e[0].append(material[i])
        material_e[1].append(material[i])
        for i in range(len(keyloc[2])):
            if element[i][0] == 'C3D8R':
                material_e[1].append([element[i][1][0][0],C3D8_pier])
            elif element[i][0] == 'S4R':
                material_e[0].append([element[i][1][0][0],element[i][1][0][-1]])
#这里是三个append是什么
material_e[0].append([500,200,1])    
material_e[2].append([2,0.1])
material_e[4].append([0.25])
#material_e[i][0]  单元的种类名称 如 'Floor'  
#material_e[i][1][0]  材料的种类 如 'CONCRETE'   
#material_e[i][1][1]  材料的密度
#material_e[i][1][2]  材料的杨氏模量和泊松比          
print('material_e',material_e)

#-----------输出ID信息 扩充ID到6维-------------
#考虑固定边界和不同单元的自由度限制
ID_ = [[1,1,0,0,0,1],[0,1,0,1,0,1],[0,0,0,1,1,1]]
# ID_[0]对应壳单元有五个自由度 所以将node_b[3],node_b[4]设置为0  对应ID_[0,0,1]
#需要用的是3,4,5 对于板单元而言 所以ID_ [1,1,0,0,0,1] 
#而如果是板单元有三个自由度也无所谓
# ID_[2]对应梁单元 自由度有6 所以将node[]均设置为0 对应ID_ [0,0,0]
#其实只需要对梁单元和板单元的自由度的ID矩阵做一个调整
#print('element[i][1][0]',element[i][0])

for i in range(len(keyloc[2])):
        for j in range(len(element[i][1][0])):
            e_num = element[i][1][0][j]   #单个单元的序号
            e_node = element[i][1][1][j]  #单个单元的节点的全局坐标
            for s in range(len(material_e)):
                if material_e[s][2][0] <= e_num <= material_e[s][2][1]:
                    for k in e_node:
                        #if node_b[k-1][5]==0:      #没有在固定边界上
                            if s==0:
                                if node_b[k-1][5]==0:
                                    node_b[k-1][0]=1 
                                    node_b[k-1][1]=1
                                    #node_b[k-1][4]=0# 板单元
                            elif s==2:
                                if node_b[k-1][5]==0:
                                    node_b[k-1][1]=1
                                    node_b[k-1][3]=1
                                    #node_b[k-1][4]=0 #梁单元
                            else:
                                if node_b[k-1][5]==0:
                                    node_b[k-1][3]=1
                                    node_b[k-1][4]=1
                                    #node_b[k-1][2]=0#其他
                        #else:        
                        #for p in range(3):
for i in range(len(node[0])):
        node_b[i][5]=1   #固定全是1
                        #    node_b[k-1][p+3]=node_b[k-1][p+3]*ID_[s][p]
for k in range(len(node_b)):                        
    if node_b[k]==[0,1,0,1,0,1]:
        print('connect',node_b[k])
                                          
#-----------输出PART信息到文件-------------   
location_write = path[2:].replace('\\','/') + '/' + 'Bridge' + '.dat' 
inp = open(location_write, 'w')

#输入名称
inp.write('The inp file for stappp' + '\n')

#输入节点号
inp.write('%10d'*4%(len(node[0]),len(keyloc[2]),1,modecontrol) + '\n')


#输入节点信息
for i in range(len(node[0])):
    inp.write('%10d'%(node[0][i]) + '%2d'*6%(node_b[i][0],node_b[i][1],node_b[i][2],node_b[i][3],node_b[i][4],node_b[i][5]))
    inp.write('%10.3f'*3%(node[1][i][0],node[1][i][1],node[1][i][2]))
    inp.write('\n')

#输入荷载信息
inp.write('%5d'*2%(1,0)+'\n') #密度在单元处给出
  

#输入单元信息
for i in range(len(keyloc[2])):

    if element[i][0] == 'C3D8R':   # 8H -4
        inp.write('%10d'*3%(4,len(element[i][1][0]),2) + '\n')
        inp.write('%10d'%1 + '%10.3e'%(25e9) + '%10.3f'%(0.3)+'%10.1f'%(2320.0) + '\n')
        inp.write('%10d'%2 + '%10.3e'%(60e9) + '%10.3f'%(0.27) +'%10.1f'%(2770.0)+ '\n')
        for j in range(len(element[i][1][0])):
            el_num = j+1
            el_node = element[i][1][1][j]
            inp.write('%10d'*9%(el_num,el_node[0],el_node[1],el_node[2],el_node[3],el_node[4],el_node[5],el_node[6],el_node[7]))
            if el_num < C3D8_pier :
                inp.write('%2d'%(1) + '\n')
            else:
                inp.write('%2d'%(2) + '\n')
    elif element[i][0] == 'S4R':  # apply plate - 6
        inp.write('%10d'*3%(6,len(element[i][1][0]),1) + '\n')
        inp.write('%10d'%1 + '%10.3e'%(25e9) + '%10.3f'%(0.3)+ '%10.3f'%(1)+'%10.1f'%(2320.0) + '\n')
        for j in range(len(element[i][1][0])):
            el_num = j+1
            el_node = element[i][1][1][j]
            inp.write('%10d'*5%(el_num,el_node[0],el_node[1],el_node[2],el_node[3]))
            inp.write('%2d'%(1)  + '\n')
    elif element[i][0] == 'B31':   # Beam -5   # Data Checked
        inp.write('%10d'*3%(5,len(element[i][1][0]),1) + '\n')
        # G=E/(2*(1+u)
        inp.write('%10d'%1 + '%10.3e'%(70e9)  + '%10.3f'%(0.346)+ '%10.1e'%(2710)+'%10.3f'%(0.76))
        inp.write('%10.3e'%(0.45853333) + '%10.3e'%(0.45853333) + '%10.3e'%(0.45853333)+'\n')
        #inp.write('%10.3e'%(0.9170666) + '%10.3e'%(0.9170666) + '%10.3e'%(0.9170666) + '\n')
        for j in range(len(element[i][1][0])):
            el_num = j+1
            el_node = element[i][1][1][j]
            inp.write('%10d'*3%(el_num,el_node[0],el_node[1]))
            inp.write('%2d'%(1) + '\n')
    elif element[i][0] == 'T3D2':  #Bar -1
        inp.write('%10d'*3%(1,len(element[i][1][0]),1) + '\n')
        inp.write('%10d'%1 + '%10.3e'%(117e9) + '%10.3f'%(0.25)+ '%10.1f'%(7860) + '\n')
        for j in range(len(element[i][1][0])):
            el_num = j+1
            el_node = element[i][1][1][j]
            inp.write('%10d'*3%(el_num,el_node[0],el_node[1]))
            inp.write('%2d'%(1) + '\n')
    
inp.write('End of INPUT file')            
inp.close()    
    
    
    
    
    
    
    
    
