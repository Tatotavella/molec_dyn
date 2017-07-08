import numpy as np
import sys
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    root_path = sys.argv[1]
else:
    root_path = './'

files_full_path = root_path + '/trajectories_*_data.csv_lindemann.npy'
files_list = glob.glob(files_full_path)

temp=[]
#extraigo las temperaturas
for i in range(0,len(files_list)):
 temp.append(float(re.findall(r"\d+\.\d+", files_list[i])[0]))

#ordeno
temp=sorted(temp)
files_list=sorted(files_list)

melt=[]
for i in range(0,len(temp)):
 data=np.load(files_list[i])
 melt.append(data[1][len(data[1])-1])
 
plt.plot(temp,melt,'-')
plt.plot(temp,melt,'bo')
plt.show()
