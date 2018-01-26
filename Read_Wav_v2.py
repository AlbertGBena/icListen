
import numpy as np
import struct
import io
import matplotlib.pyplot as plt

from pylab import *
from decimal import *
import datetime



down63=56.3
up63=70.8
down125=112
up125=141
down2k=1778
up2k=2239


filename='D:/SARTI/projectes europeus/Hidrofon/Hidrofono gallego/Rebut20180119/Datalogger_ON-CTD_OFF-ADCP_OFF/SBW1369_20180117_092200.wav'
##ifile=wave.open(filename,'rb')
Dat=[]
Dat2=[]
##Long=ifile.getnframes()
##
##Para=ifile.getparams()
##print(Long)
##print(Para)
##ifile.close()
File=io.open(filename,'rb')
bufHeader=File.read(8) #aquí he ficat els valors petits per veure després el total 
##print('complet',struct.unpack('2038s',bufHeader[4:2042]))
print('type of file',bufHeader[0:4])
print('chunk size',struct.unpack('I',bufHeader[4:8]))
Chunksize=struct.unpack('I',bufHeader[4:8])
MidaTotal=int(Chunksize[0])
File.close() #tanco arxiu i el torno a obrir amb el nombre exacte de dades


File=io.open(filename,'rb')
bufHeader=[]
bufHeader=File.read(Chunksize[0])
print('Format chunk',struct.unpack('4s',bufHeader[8:12]))

IDFile=struct.unpack('23s',bufHeader[112:135])
print('Filename',IDFile[0])
Idfile=IDFile[0].split('_')
Date=Idfile[2]+Idfile[1]
temps1=datetime.datetime.strptime(Date,'%H%M%S%Y%m%d')

Vpp2=struct.unpack('8s',bufHeader[148:156])
Vpp=int(Decimal(Vpp2[0]))
print('tensio max peak  ',Vpp)

S2=struct.unpack('6s',bufHeader[162:168])
S=int(Decimal(S2[0]))
print('sensibilitat ',S)

print('fmt  ',struct.unpack('4s',bufHeader[260:264]))
print('chunk size  ',struct.unpack('I',bufHeader[264:268]))
print(' compression code ',struct.unpack('h',bufHeader[268:270]))
print(' number of channels ',struct.unpack('h',bufHeader[270:272]))

Fs=struct.unpack('I',bufHeader[272:276])
fs=int(Fs[0])
print(' sample rate ',fs)

print(' bytes per second ',struct.unpack('I',bufHeader[276:280]))
print(' bloc alignement or bytes per sample ',struct.unpack('h',bufHeader[280:282])) #sempre que nomes hi hagi 1 canal

BRate2=struct.unpack('h',bufHeader[282:284])
BRate=int(BRate2[0])
print(' bits per sample ',BRate)

LL=struct.unpack('I',bufHeader[288:292])
print(' chunk size ',LL)


vector=[]
vector2=[]
P2=[]
P=[]
Rawdata=[]

Rawdata=bufHeader[292:int(LL[0])+292]

if BRate==24: #pot ser 16 bits or 24 bits
    DD=[]
    L24=len(Rawdata)/3
    t6=datetime.datetime.now()
    DD=np.asarray(list(''.join(l + '\x00' * (n % 3 == 2) for n, l in enumerate(Rawdata)))) #t6riga molt s'ha de optimitzar
    t7=datetime.datetime.now()
    t8=t7-t6
    print('tenmos del join',divmod(t8.total_seconds(), 60))
    vector2=np.ndarray((L24,),dtype='<i4',buffer=DD)
    vector2[vector2>8388607]-=16777216



else:
    vector2=np.ndarray((len(Rawdata)/2,),dtype='<i2',buffer=Rawdata)
    

Offset=int(round(np.mean(vector2),0))
Values2=np.subtract(vector2,Offset)

print('offset en comptes',Offset)

Values=np.array(Values2)
#convertim les contes en pressions en dB
##


Convert=20*np.log10(Vpp)-20*(BRate-1)*np.log10(2)
##print(Convert)
##agafem mostres cada 20 segons:

for j in range(1+int(round(len(vector2)/(20*fs),0))): #faig troços èr saber quants trams de 20 segons hi ha
    Data=[]
    P63=[]
    P125=[]
    P2k=[]
    Data2=[]
    Data3=[]
    temps=temps1+datetime.timedelta(seconds=20*j)
    
    if BRate==24 and j==int(round(len(vector2)/(20*fs),0)):
##        print('dintre')
        Data=vector2[fs*20*j:((j+1)*fs*20)-4] #s'ha de ficar perque si no es penja
    else:
##        print('normal')
        Data=vector2[fs*20*j:(j+1)*fs*20]
            
    Offset=int(round(np.mean(Data),0))
    Data2=np.subtract(Data,Offset)
    Data3=np.array(Data2)
    P=20*np.log10(np.abs(Data3))+Convert+168 #falta inclore el terme negatiu de la ganancia!!!!
    P2=np.power(10,P/20)
    PArms=np.power(P2,2)
    PFFT=2*np.fft.rfft(P2,norm='ortho')#es multiplica per 2 per agafar la part la energia en la part de n/2+1 fins el final
    Frq=np.fft.rfftfreq(len(Data),d=1./fs)


    Ind63=np.argwhere((Frq>=down63) & (Frq<=up63))#np.where(np.logical_and(Frq>=down63, Frq<=up63))
    Ind125=np.argwhere((Frq>=down125) & (Frq<=up125))#np.where(np.logical_and(Frq>=down63, Frq<=up63))
    Ind2k=np.argwhere((Frq>=down2k) & (Frq<=up2k))
     
    for j in range(len(Ind63)):P63.append(PFFT[Ind63[j]])
    for k in range(len(Ind125)):P125.append(PFFT[Ind125[k]])
    for m in range(len(Ind2k)):P2k.append(PFFT[Ind2k[m]])


    Prms=(np.sum(PArms)/(np.size(PArms))) # dividim LA SUMA DUACRADTICA DE TOTS ELS ELEMENTS pel nombre d'elemts
        
    SPL=round(10*np.log10(Prms),2)#i fem la seva arrel cuadrada per això es 10 en comptes de 20
 
        
    SPL63=round(10*np.log10(np.sum(np.power(np.abs(P63),2))/np.size(PFFT)),2)
    SPL125=round(10*np.log10(np.sum(np.power(np.abs(P125),2))/np.size(PFFT)),2)
    SPL2k=round(10*np.log10(np.sum(np.power(np.abs(P2k),2))/np.size(PFFT)),2)
    print(temps,'SPL in 20 seconds',SPL,'SPL in 20 seconds at 63 Hz',SPL63,'SPL in 20 seconds at 125Hz',SPL125,'SPL in 20 seconds at 2kHz',SPL2k)
##    plt.subplot(111)
##    plt.plot(Frq,np.abs(PFFT))
##
##    
##    
##    plt.show()
##    

    
                       


File.close()
