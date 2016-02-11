from rsf.proj import *

import os
env = Environment(ENV = os.environ,CC='nvcc')
env.Prepend(CPPPATH=['/home/harpreet/madagascar-1.7/include'],
			LIBPATH = ['/home/harpreet/madagascar-1.7/lib'] ,
			LIBS = ['m','rsf'] )
ldflags = '-I'+env['CPPPATH'][0]
exe=env.Program('gpumod.exe','gpumod.cu',LINKFLAGS=ldflags)
 
#Alternative ( forced )        
#exe=env.Command('gpumod.exe' , 'gpumod.cu' , 
#" nvcc $SOURCE -I/home/harpreet/madagascar-1.7/include -L/home/harpreet/madagascar-1.7/lib -lrsf -lm -o $TARGET")
#--------

Flow('vel',None,
'''
spike mag=1800,1800,2200,2200 nsp=4 k1=1,236,236,266 l1=235,265,265,500 k2=1,1,251,1 l2=500,250,500,500 o1=0 o2=0
d1=5 d2=5 n1=500 n2=500 label1=x1 unit1=m label2=x2 unit2=m 
'''
	)
Plot('vel','grey gainpanel=all color=j')
Flow('shotsa checka','vel %s' %exe[0],
	'''
	./${SOURCES[1]} check=${TARGETS[1]} csdgather=n fm=15 amp=1 dt=0.0005 ns=1 ng=500 nt=4000
	sxbeg=200 szbeg=0 jsx=0 jsz=0 gxbeg=0 gzbeg=0 jgx=1 jgz=0 chk=y kt=100
	
	''')
Plot('shotsa','grey gainpanel=all ')
Result('shotsa','grey gainpanel=all ')

Flow('velsm','vel','smooth rect2=25 repeat=10')
Plot('velsm','grey gainpanel=all color=j')
Flow('shotsb checkb','velsm %s' %exe[0],
	'''
	./${SOURCES[1]} check=${TARGETS[1]} csdgather=n fm=15 amp=1 dt=0.0005 ns=1 ng=500 nt=4000
	sxbeg=200 szbeg=0 jsx=0 jsz=0 gxbeg=0 gzbeg=0 jgx=1 jgz=0 chk=y kt=100
	
	''')
Result('velcomp','vel velsm','SideBySideAniso')
Plot('shotsb','grey gainpanel=all ')
Flow('diff','shotsa shotsb','add scale=1,-1 ${SOURCES[1]}' )
Plot('diff','sfgrey gainpanel=all')
Result('datcomp','shotsa shotsb diff','SideBySideAniso')

End()


End()
