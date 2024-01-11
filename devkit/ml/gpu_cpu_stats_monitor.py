import subprocess
from time import sleep
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

import pandas as pd


def prep(x):
    if isinstance(x,str) and '%' in x:
        return float(x.split()[0])
    elif isinstance(x,float):
        print(f'it is a float: {x}')
    return None


class NvidiaSmiLog:

    def __init__(self,n_secs=1) -> None:
        self.df = None
        self.proc = None
        self.n_secs = n_secs
        

    def start(self):
        cmd = f"nvidia-smi --query-gpu=timestamp,name,pci.bus_id,temperature.gpu,utilization.gpu,utilization.memory --format=csv -l {self.n_secs}".split()
        self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    def stop(self):
        self.proc.kill()
        output = self.proc.stdout.read()
        print(type(output),repr(output))
        csv = StringIO(output.decode())
        df = pd.read_csv(csv,sep=',')   
        df = df.rename(columns={
            ' utilization.gpu [%]':'utilization.gpu',
            ' utilization.memory [%]':'utilization.memory',
            ' temperature.gpu':'temperature.gpu',
            ' name':'name'
        })
        df['utilization.gpu'] = df['utilization.gpu'].map(prep)
        df['utilization.memory'] = df['utilization.memory'].map(prep)
        df['timestamp'] = pd.to_datetime(df['timestamp'],format= '%Y/%m/%d %H:%M:%S.%f')
        self.df = df

    def plot(self):
        dfm = pd.melt(self.df,id_vars=['timestamp','name'],value_vars=['utilization.gpu','utilization.memory'])
        sns.lineplot(x="timestamp", y="value",style='variable',hue="name",data=dfm)
        plt.ylabel('Utilization %')
        plt.xticks(rotation=90)
        ax2 = plt.twinx()
        sns.lineplot(data=self.df,x='timestamp',y='temperature.gpu', color="r", ax=ax2)
        plt.ylabel('Temp C')
        plt.xticks(rotation=90)
        plt.show()

        
"""
cmd = "nvidia-smi --query-gpu=timestamp,name,pci.bus_id,temperature.gpu,utilization.gpu,utilization.memory --format=csv -l 1".split()
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
sleep(2)
proc.kill()
output = proc.stdout.read()
print(type(output),repr(output))
TESTDATA = StringIO(output.decode())
df = pd.read_csv(TESTDATA,sep=',')
print(df)
"""
#print('-----------------')
#print(output)
