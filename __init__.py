import matplotlib.pyplot as pl
from numpy import array,arange,fromfile,empty,diff,c_
import gzip
import os
import gc

import time
from glob import glob
import re

import scipy
from scipy import interpolate
import PIL.Image
import PIL.ImageDraw
from cernlogdb import *

def makepalette(cm):
  pal=[ [r,g,b] for r,g,b,a in cm(range(256),bytes=True)]
  return reduce(lambda x,y: x+y,pal)

pal_copper=makepalette(pl.cm.copper)
pal_jet=makepalette(pl.cm.jet)

def plot2d(v,xlim=None,ylim=None,**argsv):
  if xlim is None:
    xlim=[0,len(v[0])-1]
  if ylim is None:
    ylim=[0,len(v)-1]
  p=pl.imshow(v,aspect='auto',origin='bottom',extent=xlim+ylim, **argsv)
#  pl.colorbar()
  return p

def plraw(t,v,name='raw',vmin=-160,vmax=-40):
  pl.clf()
  p=plot2d(v,xlim=[0,0.5],ylim=[0,t[-1]-t[0]],cmap=pl.cm.jet,vmin=vmin,vmax=vmax)
  #pl.yticks(  arange(0,3600,300), map(str,range(0,60,5) ))
  pl.colorbar()
  pl.title(name)
  return p

def pldetect(t,v,name='detect',vmax=2):
  v1=diff(v,axis=0)
  vmax=int(diff(t).mean())/10.*vmax
  pl.clf()
  p=plot2d(v1,vmin=-vmax,vmax=vmax,xlim=[0,0.5],ylim=[0,t[-1]-t[0]],cmap=pl.cm.copper)
  #pl.yticks(  arange(0,3600,300), map(str,range(0,60,5) ))
  pl.colorbar()
  pl.title(name)
  return p

def spraw(v,name='raw',pal=pal_jet,vmin=-160,vmax=-40):
  im=scipy.misc.toimage(v,cmin=-vmax,cmax=vmax)
  flip_h=1
  im=im.transpose(flip_h)
  if pal is not None:
    im.putpalette(pal)
  return im

def spdetect(v,name='detect',vmax=2,pal=pal_copper):
  v1=diff(v,axis=0)
  im=scipy.misc.toimage(v1,cmin=-vmax,cmax=vmax)
  flip_h=1
  im=im.transpose(flip_h)
  if pal is not None:
    im.putpalette(pal)
  return im

def mysavefig(fn):
  print fn
  pl.savefig(fn)
  


def parsetime(s):
  d=s.split('.')
  if len(d)==1:
    d=d[0];sec=0
  else:
    d=d[0];sec=float('.'+d[1])
  t=time.mktime(time.strptime(d,'%Y-%m-%d %H:%M:%S'))
  return t+sec

def dumptime(t):
  ti=int(t)
  sec=round((t-ti),3)*1000
  tt=time.localtime(ti)
  return time.strftime('%Y-%m-%d %H:%M:%S.%%d',tt)%sec

def dumptime2(t):
  ti=int(t)
  sec=round((t-ti),3)*1000
  tt=time.localtime(ti)
  return time.strftime('%Y%m%d%H%M%S',tt)

def get_chunk(vs,t1,t2,sa=None,sf=None,exe='cern-ldb'):
  """Usage get_chunk("SPS.BCTDC.31832:INT","2010-06-10 00:00:00","2010-06-10 23:59:59")
  sf  one of _1_SEC,_10_SECS,_1_MIN,_5_MIN,_10_MIN,_15_MIN,_1_HOUR,_1_DAY
  sa  one of AVG,MIN,MAX,REPEAT,INTERPOLATE,SUM,COUNT
  fmt one of CSV,XLS,TSV,MATHEMATICA
  """
  if type(t1) in [float,int]:
    t1=dumptime(t1)
  if type(t2) in [float,int]:
    t2=dumptime(t2)
  cmd='%s -vs "%s" -t1 "%s" -t2 "%s"' %(exe,vs,t1,t2)
  if sa is not None and sf is not None:
    cmd+=' -sa "%s" -sf "%s"' %(sa,sf)
  print cmd
  fh=os.popen(cmd)
  data={}
  dataon=False
  log=[]
  for l in fh.readlines():
    if dataon and l=='\n':
      dataon=False
    elif l.startswith('VARIABLE'):
      vname=l.split()[1]
      t,v=[],[]
      data[vname]=(t,v)
      dataon=False
    elif dataon:
      ll=l.strip().split(',')
      t.append(parsetime(ll[0]))
      v.append(map(float,ll[1:]))
    elif l.startswith('Timestamp'):
      dataon=True
    else:
      log.append(l)
  log=''.join(log)
  print "Found %s entries" % ','.join([ str(len(v)) for t,v in data.values()])
  return data,log

def get_data(vs,t1,t2,step,sa=None,sf=None,exe='cern-ldb'):
  if type(t1) is str:
    t1=round(parsetime(t1),0)
  if type(t2) is str:
    t2=round(parsetime(t2),0)
  datal=[]
  logl=[]
  for t in range(int(t1),int(t2),step):
    data,log=get_chunk(vs,t,t+step-0.001,sa=sa,sf=sf,exe=exe)
    datal.append( data)
    logl.append(log)
  datas={}
  for vs in data.keys():
    t,v=[],[]
    datas[vs]=(t,v)
    for data in datal:
      t.extend(data[vs][0])
      v.extend(data[vs][1])
  logs=''.join(logl)
  return datas,logs

def mk1minfn(vs,t):
  tt=time.localtime(int(t))
  return time.strftime('data1min/%%s_%Y%m%d.bin',tt)%vs

def mk10secfn(vs,t):
  tt=time.localtime(int(t))
  return time.strftime('data10sec/%%s_%Y%m%d%H.bin',tt)%vs

def mk1secfn(vs,t):
  tt=time.localtime(int(t))
  return time.strftime('data1sec/%%s_%Y%m%d%H.bin',tt)%vs

def interpdata(t,v,nt=arange(0,3600,4),havg=4):
  size=(len(nt),v.shape[1]/havg)
  img=empty(size)
  nv=v[:,0::havg]
  for i in range(1,havg):
    nv+=v[:,i::havg]
  nv/=havg
  for i in range(size[1]):
    ifun=interpolate.interp1d(t,nv[:,i],fill_value=0,bounds_error=False)
    img[:,i]=ifun(nt)
  return img

def loaddata(fn):
  vv=fromfile(fn)
#  print fn,len(vv),len(vv)%4097
  vv=vv.reshape(len(vv)/4097,4097)
  t=vv[:,0]
  v=vv[:,1:]
  return t,v

def loaddata_double(fn):
  vv=fromfile(fn)
#  print fn,len(vv),len(vv)%8193
  vv=vv.reshape(len(vv)/8193,8193)
  t=vv[:,0]
  v=vv[:,1:]
  return t,v

def loaddata_half(fn):
  vv=fromfile(fn)
#  print fn,len(vv),len(vv)%2049
  vv=vv.reshape(len(vv)/2049,2049)
  t=vv[:,0]
  v=vv[:,1:]
  return t,v


def today_hour(t=None):
  t=time.localtime(t)
  return int(parsetime(time.strftime('%Y-%m-%d %H:00:00',t)))

def today_day(t=None):
  t=time.localtime(t)
  return int(parsetime(time.strftime('%Y-%m-%d 00:00:00',t)))

b1h,b1v,b2h,b2v="LHC.BQBBQ.UA47.FFT1_B1:FFT_DATA_H LHC.BQBBQ.UA47.FFT1_B1:FFT_DATA_V LHC.BQBBQ.UA43.FFT1_B2:FFT_DATA_H LHC.BQBBQ.UA43.FFT1_B2:FFT_DATA_V".split()

