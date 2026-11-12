"""Small standard-library desktop GUI for quick and remote ANSYS contest runs."""
from __future__ import annotations
import json
import queue
import subprocess
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

ROOT=Path(__file__).resolve().parent

class PizzaApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Pizza Box Contest - Quick FEM / ANSYS via SSH')
        self.geometry('970x710')
        self.q=queue.Queue();self.running=False
        self.step=tk.StringVar();self.dxf=tk.StringVar();self.job=tk.StringVar(value=str(ROOT/'results'/'team_01'))
        self.fast_mesh=tk.StringVar(value='3.0');self.ansys_mesh=tk.StringVar(value='2.0')
        self.dt=tk.StringVar(value='10');self.sensors=tk.StringVar(value='64')
        self.host=tk.StringVar(value='');self.nproc=tk.StringVar(value='4');self.mesh_mode=tk.StringVar(value='voxel')
        self.ansysexe=tk.StringVar()
        form=ttk.Frame(self,padding=15);form.pack(fill='x')
        form.columnconfigure(1,weight=1)
        self.path(form,0,'Folded STEP (.step, .stp)',self.step,('STEP CAD','*.step *.stp'))
        self.path(form,1,'Flat DXF (.dxf; CUT/FOLD)',self.dxf,('DXF','*.dxf'))
        self.path(form,2,'Result / job directory',self.job,None)
        row=ttk.Frame(form);row.grid(row=3,column=0,columnspan=3,sticky='ew',pady=(12,5))
        for label,var,width in [('Quick mesh (mm)',self.fast_mesh,7),('ANSYS mesh (mm)',self.ansys_mesh,7),
                                ('Time step (s)',self.dt,7),('Wall samples',self.sensors,7)]:
            ttk.Label(row,text=label).pack(side='left',padx=(0,4));ttk.Entry(row,textvariable=var,width=width).pack(side='left',padx=(0,12))
        ttk.Label(row,text='ANSYS mesh type').pack(side='left',padx=(0,5))
        ttk.Combobox(row,textvariable=self.mesh_mode,values=('voxel','step'),state='readonly',width=9).pack(side='left')
        row=ttk.Frame(form);row.grid(row=4,column=0,columnspan=3,sticky='ew',pady=(5,10))
        ttk.Label(row,text='SSH user@Tailscale-IP').pack(side='left',padx=(0,6))
        ttk.Entry(row,textvariable=self.host,width=27).pack(side='left',padx=(0,12))
        ttk.Label(row,text='ANSYS CPU count').pack(side='left',padx=(0,5))
        ttk.Entry(row,textvariable=self.nproc,width=5).pack(side='left',padx=(0,12))
        ttk.Label(row,text='ANSYS exe (optional)').pack(side='left',padx=(0,5))
        ttk.Entry(row,textvariable=self.ansysexe,width=40).pack(side='left',fill='x',expand=True)
        row=ttk.Frame(form);row.grid(row=5,column=0,columnspan=3,sticky='ew')
        self.buttons=[]
        for label,fn in [('1. Quick FEM run',self.quick),('2. Prepare ANSYS',self.prepare),
                         ('3. Remote ANSYS via SSH',self.remote),('4. Score from ANSYS CSV',self.finalize),
                         ('5. Rank several teams',self.rank)]:
            b=ttk.Button(row,text=label,command=fn);b.pack(side='left',padx=4);self.buttons.append(b)
        ttk.Separator(self,orient='horizontal').pack(fill='x',padx=12,pady=5)
        status=ttk.Frame(self,padding=(15,0,15,6));status.pack(fill='x')
        self.message=tk.StringVar(value='Ready. Select a STEP and a DXF for the first two modes.')
        ttk.Label(status,textvariable=self.message).pack(side='left')
        txt=ttk.Frame(self,padding=(14,4,14,14));txt.pack(fill='both',expand=True)
        self.log=tk.Text(txt,wrap='word',font=('Menlo',11),height=20)
        scrollbar=ttk.Scrollbar(txt,orient='vertical',command=self.log.yview)
        self.log.configure(yscrollcommand=scrollbar.set)
        self.log.pack(side='left',fill='both',expand=True);scrollbar.pack(side='right',fill='y')
        self.after(150,self.tick)
    def path(self,parent,n,label,var,filetype):
        ttk.Label(parent,text=label).grid(row=n,column=0,sticky='w',pady=5)
        ttk.Entry(parent,textvariable=var).grid(row=n,column=1,sticky='ew',padx=8,pady=5)
        def browse():
            if filetype:ans=filedialog.askopenfilename(filetypes=[filetype,('All files','*.*')]);
            else:ans=filedialog.askdirectory()
            if ans:var.set(ans)
        ttk.Button(parent,text='Browse',command=browse).grid(row=n,column=2,pady=5)
    def validate_paths(self):
        if not Path(self.step.get()).is_file() or not Path(self.dxf.get()).is_file():
            messagebox.showerror('Missing input','Choose existing folded STEP and unfolded DXF files.');return False
        return True
    def execute(self,args):
        if self.running:
            messagebox.showwarning('Busy','One run at a time.');return
        self.running=True;self.message.set('Running: '+args[1])
        for b in self.buttons:b.configure(state='disabled')
        self.log.insert('end','\n$ '+subprocess.list2cmdline(args)+'\n');self.log.see('end')
        def worker():
            try:
                p=subprocess.Popen(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                   text=True,bufsize=1)
                for line in p.stdout:self.q.put(('text',line))
                code=p.wait();self.q.put(('done',code))
            except Exception as e:self.q.put(('text',f'Error: {e}\n'));self.q.put(('done',-1))
        threading.Thread(target=worker,daemon=True).start()
    def quick(self):
        if not self.validate_paths():return
        self.execute([sys.executable,str(ROOT/'simulator.py'),'--step',self.step.get(),'--dxf',self.dxf.get(),
                      '--out',self.job.get(),'--voxel-mm',self.fast_mesh.get(),'--dt',self.dt.get()])
    def prepare(self):
        if not self.validate_paths():return
        self.execute([sys.executable,str(ROOT/'ansys_bridge.py'),'prepare','--step',self.step.get(),
                      '--dxf',self.dxf.get(),'--out',self.job.get(),'--mesh-mode',self.mesh_mode.get(),
                      '--mesh-mm',self.ansys_mesh.get(),'--dt',self.dt.get(),'--sensors',self.sensors.get()])
    def remote(self):
        if not self.host.get().strip():
            messagebox.showerror('Missing SSH host','Provide username@Tailscale IP or SSH hostname.');return
        args=[sys.executable,str(ROOT/'ansys_bridge.py'),'ssh','--job',self.job.get(),
              '--host',self.host.get(),'--nproc',self.nproc.get()]
        if self.ansysexe.get().strip():args+=['--ansys-exe',self.ansysexe.get()]
        self.execute(args)
    def finalize(self):
        self.execute([sys.executable,str(ROOT/'ansys_bridge.py'),'finalize','--job',self.job.get()])
    def rank(self):
        paths=filedialog.askopenfilenames(title='Choose multiple result.json files',filetypes=[('JSON','*.json')])
        if len(paths)<2:
            messagebox.showinfo('Ranking','Choose at least two result.json files');return
        dst=filedialog.askdirectory(title='Choose ranking output directory')
        if not dst:return
        self.execute([sys.executable,str(ROOT/'ansys_bridge.py'),'rank',*paths,'--out',dst])
    def tick(self):
        try:
            while True:
                what,data=self.q.get_nowait()
                if what=='text':self.log.insert('end',data);self.log.see('end')
                elif what=='done':
                    self.running=False
                    self.message.set('Finished' if data==0 else 'Failed - inspect log')
                    for b in self.buttons:b.configure(state='normal')
                    if data==0:
                        result=Path(self.job.get())/'result.json'
                        if result.is_file():
                            try:
                                metrics=json.loads(result.read_text())['official_metrics']
                                self.log.insert('end','\nMETRICS: '+json.dumps(metrics,ensure_ascii=False,indent=2)+'\n')
                            except Exception:pass
        except queue.Empty:pass
        self.after(150,self.tick)

if __name__=='__main__':PizzaApp().mainloop()
