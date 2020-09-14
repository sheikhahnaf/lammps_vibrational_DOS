from numpy.fft import rfft, irfft
from sys import exit
from math import sqrt
from math import sin
from math import cos
from math import pi
from math import exp



########################################  Function Definition  ###################################################################

# The correlation theorem says that multiplying the Fourier transform of
# one function by the complex conjugate of the Fourier transform of the other gives the
# Fourier transform of their correlation.

def correlate(a,b):
  length = len(a)
  a = rfft(a).conjugate()     #  a(t0)b(t0+t)
  b = rfft(b)                 # .conjugate() for b(t0)a(t0+t)
  c = irfft(a*b)/length
  return c

def correlate_zhxl(a,b):  # testing correlate()
  length = len(a)
  c = [0.0 for i in range(length)]
  d = [0.0 for i in range(length)]
  for i in range(length):
    for j in range(length):
      k = i + j
      if k>=length:
        k -= length
      c[j] += a[i]*b[k]/length
  return c

########################################  Function Definition  ###################################################################

ofile1 = './VACF.normalized.dat'
ofile2 = './VDOS.normalized.dat'
ifile = './v1.txt'

fo1 = open(ofile1,'w')
fo2 = open(ofile2,'w')
fi  = open(ifile,'r')

m = 0
n = 0
read_time = False
read_natoms = False
while True:
  line = fi.readline()
  if len(line)==0:
    break
  if line.find('ITEM: TIMESTEP') != -1:
    read_time = True
    continue
  if read_time == True:
    n += 1
    a0 = int(line.split()[0])
    if n==1:
      step1 = a0
    if n==2:
      step2 = a0
      t = step2 - step1   #  Output interval
    nsteps = a0
    read_time = False
    continue
  if line.find('ITEM: NUMBER OF ATOMS') != -1:
    m += 1
    if m==1:
      read_natoms = True
    continue
  if read_natoms == True:
    a0 = int(line.split()[0])
    natoms = a0
    read_natoms = False
    continue
fi.close()

nmax = n - 1   # NO. of data

nmax = 10000

########################################  Parameters  ###################################################################

timestep = 0.0005      # ps

colx = 2
coly = 3
colz = 4

s   = 1          #  Sample interval, can be 1t, 2t, 3t, ... 
p   = nmax*t/s   #  Correlation Length, s*p is correlation time
d   = nmax*t     #  Time to output the correlation, d >= s*p
cut = nmax*t

scale    = timestep*s*p  #  normalized to 1

interval = int(d/s)
nmax     = int(nmax*t/s)
sample   = int(s/t)
cut      = int(cut/s)

########################################  Parameters  ###################################################################

ivacf2id  = [0 for i in range(natoms)]

x  = [[0.0 for i in range(nmax)] for j in range(natoms)]
y  = [[0.0 for i in range(nmax)] for j in range(natoms)]
z  = [[0.0 for i in range(nmax)] for j in range(natoms)]
c1 = [0.0 for i in range(nmax)]
c2 = [0.0 for i in range(nmax)]
c3 = [0.0 for i in range(nmax)]
vacf1 = [0.0 for i in range(nmax)]
vacf2 = [0.0 for i in range(nmax)]
vacf3 = [0.0 for i in range(nmax)]
vacf4 = [0.0 for i in range(nmax)]
vdos1 = [0.0j for i in range(nmax)]
vdos2 = [0.0j for i in range(nmax)]
vdos3 = [0.0j for i in range(nmax)]
vdos4 = [0.0j for i in range(nmax)]

fo1.write('correlation_time_ps vdos_x vdos_y vdos_z vdos_av\n')
fo2.write('Freq_THz vdos_x vdos_y vdos_z vdos_av\n')

fi = open(ifile,'r')
n = -1
q = -1
id = -1
read_data = False
while True:
  line = fi.readline()
  if len(line)==0:
    break
  if line.find('ITEM: TIMESTEP') != -1:
    id = -1
    read_data = False
    continue
  if line.find('ITEM: ATOMS') != -1:
    read_data = True
    n += 1
    continue
  if read_data == True:
    if n==0:
      a0 = int(line.split()[0])
      id += 1
      ivacf2id[id] = a0
      continue
    if n>0 and n%sample==0:
      id += 1
      if id==0:
        q += 1
      a0 = int(line.split()[0])
      a1 = float(line.split()[colx])
      a2 = float(line.split()[coly])
      a3 = float(line.split()[colz])
      for i in range(natoms):
        if ivacf2id[i]==a0:
          real_id = i
      x[real_id][q] = a1
      y[real_id][q] = a2
      z[real_id][q] = a3
      if id==natoms-1 and (q+1)%interval==0:
        fo1.write('%8i  '   %((q+1)*s))
        fo1.write('%8i  '   %p)
        fo1.write('\n')
        fo2.write('%8i  '   %((q+1)*s))
        fo2.write('%8i  '   %p)
        fo2.write('\n')
        print('Reading ', (q+1)*s)
        for j in range(q+1):
          vacf1[j] = 0.0
          vacf2[j] = 0.0
          vacf3[j] = 0.0
          vacf4[j] = 0.0
        for i in range(natoms):
          x1  = [0.0 for j in range(q+1)]
          y1  = [0.0 for j in range(q+1)]
          z1  = [0.0 for j in range(q+1)]
          for j in range(q+1):
            x1[j] = x[i][j]
            y1[j] = y[i][j]
            z1[j] = z[i][j]
          c1 = correlate(x1,x1)
          c2 = correlate(y1,y1)
          c3 = correlate(z1,z1)
          for j in range(q+1):
            vacf1[j] += c1[j]
            vacf2[j] += c2[j]
            vacf3[j] += c3[j]
            vacf4[j] += (c1[j]+c2[j]+c3[j])/3.0
        vacf10 = vacf1[0]
        vacf20 = vacf2[0]
        vacf30 = vacf3[0]
        vacf40 = vacf4[0]
        for j in range(q+1):
          vacf1[j] /= vacf10
          vacf2[j] /= vacf20
          vacf3[j] /= vacf30
          vacf4[j] /= vacf40
        for m in range(int(p)):
          fo1.write('%20.15f '  %((m+1)*s*timestep))
          fo1.write('%20.15f '  %vacf1[m])
          fo1.write('%20.15f '  %vacf2[m])
          fo1.write('%20.15f '  %vacf3[m])
          fo1.write('%20.15f '  %vacf4[m])
          fo1.write('\n')
        vdos1 = rfft(vacf1,q+1)
        vdos2 = rfft(vacf2,q+1)
        vdos3 = rfft(vacf3,q+1)
        vdos4 = rfft(vacf4,q+1)
        for m in range(int(p/2+1)):
          vdosx   = sqrt((vdos1[m].real)**2+(vdos1[m].imag)**2)*2.0/float(p) * scale
          vdosy   = sqrt((vdos2[m].real)**2+(vdos2[m].imag)**2)*2.0/float(p) * scale
          vdosz   = sqrt((vdos3[m].real)**2+(vdos3[m].imag)**2)*2.0/float(p) * scale
          vdos_av = sqrt((vdos4[m].real)**2+(vdos4[m].imag)**2)*2.0/float(p) * scale
          freq = float(m)/float(p)/float(s)/timestep
          fo2.write('%20.15f '  %freq)  # THz
          fo2.write('%20.15f '  %vdosx)
          fo2.write('%20.15f '  %vdosy)
          fo2.write('%20.15f '  %vdosz)
          fo2.write('%20.15f '  %vdos_av)
          fo2.write('\n')
        if q+1==cut:
          break
    continue
fi.close()
fo1.close()
fo2.close()

print('VACF and VDOS caculated OK')
