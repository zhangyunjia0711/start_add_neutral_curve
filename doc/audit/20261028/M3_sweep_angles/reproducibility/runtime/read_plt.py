"""Restricted ordered-zone v112 reader, per installed Tecplot data-format guide.

Fails explicitly on unsupported types; validates markers, extents and min/max.
"""
import struct
import numpy as np
class Reader:
    def __init__(self,b):
        self.b=b;self.p=8;assert b[:8] in [b'#!TDV112',b'#!TDV191']
        self.e='<' if b[8:12]==struct.pack('<i',1) else '>';assert self.i()==1
        self.filetype=self.i();self.title=self.s();nv=self.i();self.variables=[self.s() for _ in range(nv)];self.zones=[];self.aux={}
        while True:
            marker=self.f()
            if marker==357:break
            if marker in [799,899]:
                if marker==899:self.i()
                k=self.s();assert self.i()==0;self.aux[k]=self.s();continue
            assert marker==299,('Unsupported header marker',marker,self.p)
            z=dict(name=self.s(),parent=self.i(),strand=self.i(),time=self.d(),color=self.i(),kind=self.i());assert z['kind']==0,z
            locations=[self.i() for _ in range(nv)] if self.i() else [0]*nv
            assert not any(locations),'Cell centered not supported'
            assert self.i()==0 and self.i()==0,'Face connectivity not supported'
            z['shape']=tuple(self.i() for _ in range(3));assert min(z['shape'])>0
            z['aux']={}
            while self.i():
                k=self.s();assert self.i()==0;z['aux'][k]=self.s()
            self.zones.append(z)
        self.data_offset=self.p
    def i(self):v=struct.unpack_from(self.e+'i',self.b,self.p)[0];self.p+=4;return v
    def f(self):v=struct.unpack_from(self.e+'f',self.b,self.p)[0];self.p+=4;return v
    def d(self):v=struct.unpack_from(self.e+'d',self.b,self.p)[0];self.p+=8;return v
    def s(self):
        out=[]
        while True:
            v=self.i()
            if v==0:return ''.join(out)
            assert 0<=v<65536;out.append(chr(v))
    def data(self):
        nv=len(self.variables);self.p=self.data_offset
        for zi,z in enumerate(self.zones):
            assert self.f()==299
            types=[self.i() for _ in range(nv)];passive=[self.i() for _ in range(nv)] if self.i() else [0]*nv
            sharing=[self.i() for _ in range(nv)] if self.i() else [-1]*nv
            assert self.i()==-1
            ranges=[None if passive[v] or sharing[v]>=0 else (self.d(),self.d()) for v in range(nv)]
            arrays=[];n=int(np.prod(z['shape']))
            for v in range(nv):
                if passive[v]:a=None
                elif sharing[v]>=0:
                    assert sharing[v]<zi;a=self.zones[sharing[v]]['arrays'][v];assert len(a)==n
                else:
                    dt=np.dtype(self.e+{1:'f4',2:'f8',3:'i4',4:'i2',5:'u1'}[types[v]])
                    a=np.frombuffer(self.b,dtype=dt,count=n,offset=self.p);self.p+=n*dt.itemsize
                    finite=a[np.isfinite(a)]
                    if len(finite):
                        lo,hi=ranges[v];assert np.isclose(finite.min(),lo,rtol=2e-6,atol=1e-15) and np.isclose(finite.max(),hi,rtol=2e-6,atol=1e-15),(zi,v,ranges[v],finite.min(),finite.max())
                arrays.append(a)
            z.update(arrays=arrays,passive=passive,sharing=sharing,ranges=ranges)
        assert self.p==len(self.b),('Trailing bytes',self.p,len(self.b))
        return self
