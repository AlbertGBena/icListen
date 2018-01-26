import numpy as np
import glob
import time
import struct
import io
from decimal import *
import datetime as dt
from datetime import datetime




def diff(list1, list2):
    c = set(list1).union(set(list2))  
    d = set(list1).intersection(set(list2))
    return map(lambda s: s.strip(),sorted(list(c - d)))




down63=56.3
up63=70.8
down125=112
up125=141
down2k=1778
up2k=2239

def ReadWav(fileid,ROOT):

    filename=fileid

    Dat=[]
    Dat2=[]

    File=io.open(filename,'rb')
    bufHeader=File.read(8) 
    Chunksize=struct.unpack('I',bufHeader[4:8])
    MidaTotal=int(Chunksize[0])
    File.close() 


    File=io.open(filename,'rb')
    bufHeader=[]
    bufHeader=File.read(Chunksize[0])
    

    IDFile=struct.unpack('23s',bufHeader[112:135])
    
    Idfile=IDFile[0].split('_')
    Date=Idfile[2]+Idfile[1]
    temps1=datetime.strptime(Date,'%H%M%S%Y%m%d')

    Vpp2=struct.unpack('8s',bufHeader[148:156])
    Vpp=int(Decimal(Vpp2[0]))
    

    S2=struct.unpack('6s',bufHeader[162:168])
    S=int(Decimal(S2[0]))
    

    Fs=struct.unpack('I',bufHeader[272:276])
    fs=int(Fs[0])
    
    BRate2=struct.unpack('h',bufHeader[282:284])
    BRate=int(BRate2[0])
   

    LL=struct.unpack('I',bufHeader[288:292])
   


    vector=[]
    vector2=[]
    P2=[]
    P=[]
    Rawdata=[]

    Rawdata=bufHeader[292:int(LL[0])+292]

    if BRate==24: 
        DD=[]
        L24=len(Rawdata)/3
        t6=datetime.now()
        DD=np.asarray(list(''.join(l + '\x00' * (n % 3 == 2) for n, l in enumerate(Rawdata))))
        t7=datetime.now()
        t8=t7-t6
    
        vector2=np.ndarray((L24,),dtype='<i4',buffer=DD)
        vector2[vector2>8388607]-=16777216



    else:
        vector2=np.ndarray((len(Rawdata)/2,),dtype='<i2',buffer=Rawdata)
        

    Offset=int(round(np.mean(vector2),0))
    Values2=np.subtract(vector2,Offset)

    

    Values=np.array(Values2)
    
    

    Convert=20*np.log10(Vpp)-20*(BRate-1)*np.log10(2)
    

    for j in range(1+int(round(len(vector2)/(20*fs),0))): 
        Data=[]
        P63=[]
        P125=[]
        P2k=[]
        Data2=[]
        Data3=[]
        temps=temps1+dt.timedelta(seconds=20*j)
        
        if BRate==24 and j==int(round(len(vector2)/(20*fs),0)):
   
            Data=vector2[fs*20*j:((j+1)*fs*20)-4] 
        else:
    
            Data=vector2[fs*20*j:(j+1)*fs*20]
                
        Offset=int(round(np.mean(Data),0))
        Data2=np.subtract(Data,Offset)
        Data3=np.array(Data2)
        P=20*np.log10(np.abs(Data3))+Convert+168 
        P2=np.power(10,P/20)
        PArms=np.power(P2,2)
        PFFT=2*np.fft.rfft(P2,norm='ortho')
        Frq=np.fft.rfftfreq(len(Data),d=1./fs)


        Ind63=np.argwhere((Frq>=down63) & (Frq<=up63))
        Ind125=np.argwhere((Frq>=down125) & (Frq<=up125))
        Ind2k=np.argwhere((Frq>=down2k) & (Frq<=up2k))
         
        for j in range(len(Ind63)):P63.append(PFFT[Ind63[j]])
        for k in range(len(Ind125)):P125.append(PFFT[Ind125[k]])
        for m in range(len(Ind2k)):P2k.append(PFFT[Ind2k[m]])


        Prms=(np.sum(PArms)/(np.size(PArms))) 
            
        SPL=round(10*np.log10(Prms),2)
     
            
        SPL63=round(10*np.log10(np.sum(np.power(np.abs(P63),2))/np.size(PFFT)),2)
        SPL125=round(10*np.log10(np.sum(np.power(np.abs(P125),2))/np.size(PFFT)),2)
        SPL2k=round(10*np.log10(np.sum(np.power(np.abs(P2k),2))/np.size(PFFT)),2)
        
        temps2=temps1.strftime('%y%m%d')
        FID=open(ROOT+str(temps2)+'.csv','a+')
        FID.write(str(temps)+','+str(SPL)+','+str(SPL63)+','+str(SPL125)+','+str(SPL2k)+'\n')
        FID.close()
        


    File.close()


c=[]


fid='/media/usb0/icListen/'
Root=fid

while True:
    
    r=0
    t=datetime.now().strftime('%M') 
    if t=='30' or t=='00': 
        fid1=open(fid+'Sent.txt','a+') 

        Tot=fid1.readlines() 
        TotOrdenat=sorted(Tot)
        

        fid1.close() 
        fid1=open(fid+'Sent.txt','a')
        d=glob.glob(fid+'*.wav')
        d=sorted(d)

        n=map(lambda s: s.strip(),TotOrdenat)

        while r==0:
            
            f=diff(n,d)

            if len(f)==0:
                r=1

            else:
                n.append(f[0])
                fid1.write(f[0]+'\n')
                ReadWav(f[0],Root)
                
        n=[]
        fid1.close()

        

                
            
        

    
