import subprocess
from time import sleep
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import psutil
import threading
from threading import Thread
import tempfile
import pandas as pd
from datetime import datetime
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


class ThreadWithReturnValue(Thread):
    
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args,
                                                **self._kwargs)
    def join(self, *args):
        Thread.join(self, *args)
        return self._return

class ComputeResourceLog:

    def __init__(self,interval=1) -> None:
        self.interval = interval
        self.stop_event = threading.Event()
        self.loop_thread = ThreadWithReturnValue(target=self._log_loop)
        self.log = []
        pass

    @staticmethod
    def _get_gpu():
        cmd = f"nvidia-smi --query-gpu=timestamp,name,temperature.gpu,utilization.gpu,utilization.memory --format=csv,noheader,nounits".split()
        result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
        return result.stdout
    
    @staticmethod
    def _get_cpu():
        now = datetime.now()
        a = psutil.virtual_memory().percent
        b = psutil.virtual_memory().available * 100 / psutil.virtual_memory().total
        c = psutil.cpu_percent()
        t = now.strftime('%Y/%m/%d %H:%M:%S.%f')
        return f"{t},{a},{b},{c}\n"
    
    def _gpu_log_to_df(self,fname):
        df = pd.read_csv(fname,sep=',',names=['timestamp','name','temperature.gpu','utilization.gpu','utilization.memory'])
        df['timestamp'] = pd.to_datetime(df['timestamp'],format= '%Y/%m/%d %H:%M:%S.%f')
        return df

    def _cpu_log_to_df(self,fname):
        df = pd.read_csv(fname,sep=',',names=['timestamp','vram','vram2','cpu'])
        df['timestamp'] = pd.to_datetime(df['timestamp'],format= '%Y/%m/%d %H:%M:%S.%f')
        return df
    
    def _log_loop(self):
        temp_dir = tempfile.mkdtemp()        
        temp_file1 = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
        temp_file2 = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)

        with open(temp_file1.name, 'a') as gpu,open(temp_file2.name, 'a') as cpu:
            while not self.stop_event.is_set():
                gpu.write(self._get_gpu())
                cpu.write(self._get_cpu())
                sleep(self.interval)

        return temp_file1.name,temp_file2.name
    
    def start(self):
        self.loop_thread.start()

    def stop(self):
        self.stop_event.set()
        gpu_fname,cpu_fname = self.loop_thread.join()
        self.gpu_df = self._gpu_log_to_df(gpu_fname) 
        self.cpu_df = self._cpu_log_to_df(cpu_fname)

    def cpu_plot(self):
        plt.figure(figsize=(12,4))
        dfm = pd.melt(self.cpu_df,id_vars=['timestamp'],value_vars=['vram','vram2','cpu'])
        sns.lineplot(x="timestamp", y="value",hue='variable',data=dfm)
        plt.ylabel('Utilization %')
        plt.xticks(rotation=90)
        plt.legend(bbox_to_anchor=(1.4, 1.05))
        
        for dt,msg in self.log:
            plt.axvline(dt,color='k',linestyle='--',linewidth=0.5)
            plt.text(dt,0,msg,rotation=90)

    
    def gpu_plot(self):
        plt.figure(figsize=(12,4))

        dfm = pd.melt(self.gpu_df,id_vars=['timestamp','name'],value_vars=['utilization.gpu','utilization.memory'])
        sns.lineplot(x="timestamp", y="value",style='variable',hue="name",data=dfm)
        plt.ylabel('Utilization %')
        plt.xticks(rotation=90)
        plt.legend(bbox_to_anchor=(1.4, 1.05))
        
        for dt,msg in self.log:
            plt.axvline(dt,color='k',linestyle='--',linewidth=0.5)
            plt.text(dt,0,msg,rotation=90)
            

        ax2 = plt.twinx()
        sns.lineplot(data=self.gpu_df,x='timestamp',y='temperature.gpu', color="r", ax=ax2)
        plt.ylabel('Temp C')
        plt.xticks(rotation=90)
        plt.show() 

    def log_event(self,msg):
        now = datetime.now()
        self.log.append((now,msg))