from threading import Thread, Lock

import os
import os.path
import time
import traceback
from serial import Serial, SerialException


legs=[None,None,None,None]


WAITMOVE_DST=[50.0]

####### COMPUTE ACTUATOR ORDERS ########
def to_linear_actuator_order(l):
    """
    transform a list of elongation in an encoded msg (to send to tell_controler())
    :param l: list of 3 elongations
    :return: encoded msg
    """
    k=[(c-0.45)*1000.0 for c in l] # m to mm
    return "A%.2f#B%.2f#C%.2f#\n"%(k[2],k[0],k[1])


from subprocess import Popen,PIPE
import select
class SimuControler:
    def __init__(self,id):
        self.id=id
        self.port="/dev/simu_%d"%(id)
        self.process=Popen(['./verrin/simulator',str(id)],stdin=PIPE,stdout=PIPE)
    def isOpen(self):
        return self.process.poll()==None
    def close(self):
        self.process.kill()
        self.process.wait()
    def write(self,buffer):
        self.process.stdin.write(buffer.encode())
        self.process.stdin.flush()
    def read(self,size):    
        r=select.select([self.process.stdout.fileno()],[],[],0.01)
        if (len(r[0])>0):
            s=self.process.stdout.read(size)
            return s.decode()
        return ''
    def reset_input_buffer(self):  
        self.process.stdin.flush()
        self.process.stdout.flush()

class FakeControler:
    def __init__(self,id):
        self.id=id
        self.port="/dev/fake_%d"%(id)
        self.response=''
        self.last=time.time()
        self.la=[{'enable':False,'position':0,'order':0},{'enable':False,'position':0,'order':0},{'enable':False,'position':0,'order':0}]
    def isOpen(self):
        return True
    def close(self):
        pass
    def write(self,buffer):
        if buffer.find('?')!=-1:
            self.response+="id[%d]\r\n"%(self.id)
        else:
            d=buffer.split('#')
            for o in d:
                if len(o)>0:
                    i=ord(o[0])-ord('A')
                    if o[1:]!='_':
                        self.la[i]['enable']=True
                        self.la[i]['order']=float(o[1:])
                    else:
                        self.la[i]['enable']=False
                        #print "Fake receive: %s"%(buffer)
    def read(self,size):
        time.sleep(0.1)
        if (time.time()-self.last)>0.2:
            for l in self.la:
                if l['enable']:
                    if abs(l['order']-l['position'])<1:
                        l['enable']=False
                    else:
                        direction=-(l['position'] - l['order']) / abs(l['position']-l['order'])
                        stength=min(32,abs(l['position']-l['order']))
                        l['position'] += direction * (time.time()-self.last)*stength
            self.last=time.time()
            self.response+="S %d;%lf;%lf;0 %d;%lf;%lf;0 %d;%lf;%lf;0 \r\n"%(self.la[0]['enable'],self.la[0]['position'],self.la[0]['order'],
                                                                            self.la[1]['enable'],self.la[1]['position'],self.la[1]['order'],
                                                                            self.la[2]['enable'],self.la[2]['position'],self.la[2]['order'])
        msg=self.response[:size]
        self.response = self.response[size:]
        return msg
    def reset_input_buffer(self):
        pass


def read_id(s):
    print("read id in : ",s)
    i=s.find("id[")
    if i==-1:
        return -1
    j=s[i+len("id[")]
    return int(j)


find_mutex = Lock()

def find_controler():
    if find_mutex.acquire(False)==False:
        return
    blacklist=[] # port already up and ready
    for i in range(len(legs)):
        if legs[i]!=None and legs[i].isOpen():
            try:
                #legs[i].write('?')
                blacklist.append(legs[i].port)
            except:
                legs[i].close()
                legs[i]=None

    
    for i in range(20):
        pname="/dev/ttyUSB%d"%i
        print("try ",pname)
        if pname in blacklist:
            continue
        if os.path.exists(pname)==False:
            continue
        try:
            port=Serial(port=pname, baudrate=115200, timeout=0.2)
            print('port ',pname,' opened')
        except Exception as e:
            traceback.print_exc()
            continue
        c=0
        while port.isOpen()==False and c<10:
            time.sleep(0.5)
            c+=1
        if c==10:
            continue
        time.sleep(1.5)
        print('retrieve id...')
        port.write("???")
        time.sleep(1)
        port.write("???")
        t=0
        j=-1
        while t<10 and j==-1:
            j=read_id(port.read(8000))
        print("find id: ",j)
        if j<0 or j>3:
            continue
        if legs[j]!=None:
            legs[j].close()
        print("set port ",port.port," for leg ",j)
        legs[j]=port

    for i in range(len(legs)):
        if legs[i]==None:
            print("ADD FAKE CONTROLER!")
            legs[i]=SimuControler(i)
    find_mutex.release()


def tell_controler(number,msg,w=0):
    """
    Write encoded 'msg' in the 'number' controler

    :param number: leg id
    :param msg: encoded msg with 'to_linear_activator_order'
    :param w: waiting time after telling the controler
    :return:
    """
    print('tell controler ',number,' : ',msg)
    if legs[number]==None:
        return
    try:
        legs[number].write(msg)
        if w>0:
            wait_move(number,w)
    except Exception as e:
        traceback.print_exc()
        find_controler()
        time.sleep(0.5)
        tell_controler(number,msg,w)

def tell_controlers(V):
    """
    Send the encoded msg to all controlers from a list a 12 elongations

    :param V: list of 12 elongations
    :return:
    """
    for i in range(4):
        tell_controler(i, to_linear_actuator_order([V[3 * i] / 1000, V[3 * i + 1] / 1000, V[3 * i + 2] / 1000]))

class ControlerHandler(Thread):
    def __init__(self,id):
        Thread.__init__(self)
        self.id=id
        self.indata=''
        self.shouldStop=False
        self.isSuspended=False
        self.la=[{'enabled':False,'position':0,'target':0,'last_order':0,'time':0},
                 {'enabled':False,'position':0,'target':0,'last_order':0,'time':0},
                 {'enabled':False,'position':0,'target':0,'last_order':0,'time':0}]

    def parse_msg(self,msg):
#        if self.id==0:
#            print "parse <<:"
#            print self.indata
#            print ">>"
        if len(msg)<=0:
            return
        if msg[0]=='S':
            d=msg.split()
            if len(d)!=4:
                #print "error in data: ",d
                return
            la=0
            for v in d[1:]:
                x=v.split(';')
                if len(x)!=4:
                    print("error in data values: ",x)
                    return
                if x[0]=='0':
                    self.la[la]['enabled']=False
                else:
                    self.la[la]['enabled']=True
                self.la[la]['position'] = float(x[1])
                self.la[la]['target'] = float(x[2])
                self.la[la]['last_order'] = float(x[3])
                self.la[la]['time']=time.time()
                la+=1
            if la!=3:
                print("missing actuator in status for ",self.id)
        elif msg[0]=='D':
            print("[Debug:",self.id,"][",msg,"]")


    def print_actuators_status(self):
        print(self.id,self.la)

    def moving(self,dst):
        for c in self.la:
            if c['enabled'] and  abs(c['target']-c['position'])>=dst:
                return True
            if (time.time()-c['time'])>2:
                print("lost status for leg ",self.id)
        return False

    def stop(self):
        self.shouldStop=True
        
    def suspend(self,status):
        if status==False: # purge pending data            
            legs[self.id].reset_input_buffer()
        self.isSuspended=status
        
    def run(self):
        while self.shouldStop==False:
            try:
                if self.isSuspended:
                    time.sleep(0.5)
                    continue
                if legs[self.id]==None:
                    print("no com")
                    time.sleep(2)
                    continue
                c=legs[self.id].read(60)
                if len(c)>0:
                    full_msg=''
                    for d in c:
                        if d=='\r':
                            if len(self.indata)>0:
                                full_msg=self.indata
                                self.indata=''
                        elif d!='\n':
                            self.indata+= str(d)
                    if full_msg!='':
                        self.parse_msg(full_msg)
            except SerialException as e:
                legs[self.id]=None
                find_controler()
            except Exception as e :
                print("thread ",self.id," receive an exception:",e)
                traceback.print_exc()
        if legs[self.id]!=None:
            legs[self.id].close()
        

def wait_move(cl,timeout,dst=None):
    global WAITMOVE_DST
    if dst==None:
        dst=WAITMOVE_DST[0]
    if type(cl)==int:
        cl=(cl,)
    n=time.time()
    time.sleep(0.8) # time for order to go to physical controler and get back
    done=False
    while done==False:
        done=True
        for c in cl:
            done=done and (controlers[c].moving(dst)==False)
        print("done => ",done)
        time.sleep(0.1)
        if done==False and (time.time()-n)>timeout:
            print("wait end by timeout")
            for c in cl:
                print("leg ",c," moving:",controlers[c].moving(dst))
                if controlers[c].moving(dst):
                    for l in controlers[c].la:
                        print(l['enabled'],'/',abs(l['position']-l['target']),',', end=' ')
                print() 
                    
                #    controlers[c].print_actuators_status()
            return False
    print("wait end as there is no move",dst,": ", end=' ')
    for c in cl:
        print(c,'moving:',controlers[c].moving(dst),'[', end=' ')
        for l in controlers[c].la:
            print(l['enabled'],'/',abs(l['position']-l['target']),',', end=' ')
        print(']')
    print()
    return True
    



def upgrade_controlers(cl=None):
    if cl==None:
        cl=list(range(4))
    for c in controlers:
        c.suspend(True)
    time.sleep(0.5)
    for i in cl:
        os.system("cd verrin && CONTROLER_ID=%d ARDUINO_PORT=%s make clean all upload"%(i,legs[i].port))
    for c in controlers:
        c.suspend(False)

def stop_all_actuators():
    for i in range(4):
        tell_controler(i,"A_#B_#C_#")


print("find controlers")        
find_controler()
print("done")


print("run controler handlers")
controlers=[]
for i in range(4):
    controlers.append(ControlerHandler(i))
    controlers[-1].start()

print("done")
