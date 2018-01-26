import numpy as np
import glob
import time
##from datetime import datetime,now
import struct
import io
##import datetime as datetime
from decimal import *
import datetime as dt
from datetime import datetime


##test de proba

def diff(list1, list2):
    c = set(list1).union(set(list2))  # or c = set(list1) | set(list2)
    d = set(list1).intersection(set(list2))  # or d = set(list1) & set(list2)
    return map(lambda s: s.strip(),sorted(list(c - d)))

##funcion de leer wav


down63=56.3
up63=70.8
down125=112
up125=141
down2k=1778
up2k=2239

def ReadWav(fileid):
##    print(fileid)
    filename=fileid
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
    ##print('type of file',bufHeader[0:4])
    ##print('chunk size',struct.unpack('I',bufHeader[4:8]))
    Chunksize=struct.unpack('I',bufHeader[4:8])
    MidaTotal=int(Chunksize[0])
    File.close() #tanco arxiu i el torno a obrir amb el nombre exacte de dades


    File=io.open(filename,'rb')
    bufHeader=[]
    bufHeader=File.read(Chunksize[0])
    ##print('Format chunk',struct.unpack('4s',bufHeader[8:12]))

    IDFile=struct.unpack('23s',bufHeader[112:135])
    ##print('Filename',IDFile[0])
    Idfile=IDFile[0].split('_')
    Date=Idfile[2]+Idfile[1]
    temps1=datetime.strptime(Date,'%H%M%S%Y%m%d')

    Vpp2=struct.unpack('8s',bufHeader[148:156])
    Vpp=int(Decimal(Vpp2[0]))
    ##print('tensio max peak  ',Vpp)

    S2=struct.unpack('6s',bufHeader[162:168])
    S=int(Decimal(S2[0]))
    ##print('sensibilitat ',S)
    ##
    ##print('fmt  ',struct.unpack('4s',bufHeader[260:264]))
    ##print('chunk size  ',struct.unpack('I',bufHeader[264:268]))
    ##print(' compression code ',struct.unpack('h',bufHeader[268:270]))
    ##print(' number of channels ',struct.unpack('h',bufHeader[270:272]))

    Fs=struct.unpack('I',bufHeader[272:276])
    fs=int(Fs[0])
    ##print(' sample rate ',fs)
    ##
    ##print(' bytes per second ',struct.unpack('I',bufHeader[276:280]))
    ##print(' bloc alignement or bytes per sample ',struct.unpack('h',bufHeader[280:282])) #sempre que nomes hi hagi 1 canal

    BRate2=struct.unpack('h',bufHeader[282:284])
    BRate=int(BRate2[0])
    ##print(' bits per sample ',BRate)

    LL=struct.unpack('I',bufHeader[288:292])
    ##print(' chunk size ',LL)


    vector=[]
    vector2=[]
    P2=[]
    P=[]
    Rawdata=[]

    Rawdata=bufHeader[292:int(LL[0])+292]

    if BRate==24: #pot ser 16 bits or 24 bits
        DD=[]
        L24=len(Rawdata)/3
        t6=datetime.now()
        DD=np.asarray(list(''.join(l + '\x00' * (n % 3 == 2) for n, l in enumerate(Rawdata)))) #triga molt s'ha de optimitzar
        t7=datetime.now()
        t8=t7-t6
    ##    print('tenmos del join',divmod(t8.total_seconds(), 60))
        vector2=np.ndarray((L24,),dtype='<i4',buffer=DD)
        vector2[vector2>8388607]-=16777216



    else:
        vector2=np.ndarray((len(Rawdata)/2,),dtype='<i2',buffer=Rawdata)
        

    Offset=int(round(np.mean(vector2),0))
    Values2=np.subtract(vector2,Offset)

    ##print('offset en comptes',Offset)

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
        temps=temps1+dt.timedelta(seconds=20*j)
        
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
        Root='C://Users/Albert GB/Documents/fILES/'
        temps2=temps1.strftime('%y%m%d')
        FID=open(Root+str(temps2)+'.csv','a+')
        FID.write(str(temps)+','+str(SPL)+','+str(SPL63)+','+str(SPL125)+','+str(SPL2k)+'\n')
        FID.close()
        


    File.close()


c=[]

##bucle per identificar el nom dels arxius
fid='C://Users/Albert GB/Documents/fILES'
##d=glob.glob(fid+'/*.wav')
##print(d)

while True:
    
    r=0
    t=datetime.now().strftime('%S') #para encontrar marca de tiempo para empezar a ahacer el procesado cada 30 minutos
    if t=='30' or t=='00': #
        fid1=open(fid+'/Sent.txt','a+') #abrimos archivo de texto

        Tot=fid1.readlines() #la variable Tot contiene todos los nombres
        TotOrdenat=sorted(Tot) #ordenamos por nombre
        

        fid1.close() #cerramos archivo de texto
        fid1=open(fid+'/Sent.txt','a')
        d=glob.glob(fid+'/*.wav')
        d=sorted(d)

        n=map(lambda s: s.strip(),TotOrdenat)

        
        while r==0:
            
            f=diff(n,d)

            if len(f)==0:
                r=1

            else:
                n.append(f[0])
                fid1.write(f[0]+'\n')
                ReadWav(f[0])
                
        n=[]
        fid1.close()

        

                
            
        

    
