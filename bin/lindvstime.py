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

#ordeno
files_list=sorted(files_list)

temp=[]
#extraigo las temperaturas
for i in range(0,len(files_list)):
 temp.append(float(re.findall(r"\d+\.\d+", files_list[i])[0]))

for idx in range(0,len(files_list),2):
 filename=files_list[idx]
 data=np.load(filename)
 # grafico
 plt.plot(data[0], data[1], '-',label='%.2f' % (temp[idx]))
 plt.grid()
plt.legend(loc='lower right')
plt.grid()
plt.show()

