from serial import *
import time
import sys

leg0=Serial(port="/dev/ttyUSB0", baudrate=115200)
leg1=Serial(port="/dev/ttyUSB1", baudrate=115200)
leg2=Serial(port="/dev/ttyUSB2", baudrate=115200)
leg3=Serial(port="/dev/ttyUSB3", baudrate=115200)

c=0
while leg0.isOpen()==False or leg1.isOpen()==False or leg2.isOpen()==False or leg3.isOpen()==False:
    time.sleep(1);
    if c==10:
        sys.exit(2);
    c+=1

print("all opened")
time.sleep(5)
# print leg0.write("A100#B200#C70#")
# print leg1.write(b"A100#B200#C70#")
# print leg2.write(b"A100#B200#C70#")
# print leg3.write(b"A100#B200#C70#")

# time.sleep(3.5)

# print leg0.write("A100#B0#C120#")
# print leg1.write(b"A100#B0#C120#")
# print leg2.write(b"A100#B0#C120#")
# print leg3.write(b"A100#B0#C120#")

# time.sleep(3.5)
# print leg0.write("A100#B200#C70#")
# print leg1.write(b"A100#B200#C70#")
# print leg2.write(b"A100#B200#C70#")
# print leg3.write(b"A100#B200#C70#")

# time.sleep(3.5)

# print leg0.write("A0#B200#C70#")
# print leg1.write(b"A0#B200#C70#")
# print leg2.write(b"A0#B200#C70#")
# print leg3.write(b"A0#B200#C70#")

# time.sleep(3.5)
# print leg0.write("A200#B200#C70#")
# print leg1.write(b"A200#B200#C70#")
# print leg2.write(b"A200#B200#C70#")
# print leg3.write(b"A200#B200#C70#")

# time.sleep(3.5)

# print leg0.write("A100#B200#C70#")
# print leg1.write(b"A100#B200#C70#")
# print leg2.write(b"A100#B200#C70#")
# print leg3.write(b"A100#B200#C70#")

# time.sleep(3.5)


while True:
    m=input("?")
    if m[0]=='0':
        leg0.write(m[1:])
    elif m[0]=='1':
        leg1.write(m[1:])
    elif m[0]=='2':
        leg2.write(m[1:])
    else :
        leg3.write(m[1:])
    
    

leg0.write("A_#B_#C_#")
leg1.write("A_#B_#C_#")
leg2.write("A_#B_#C_#")
leg3.write("A_#B_#C_#")

leg0.close()
leg1.close()
leg2.close()
leg3.close()
