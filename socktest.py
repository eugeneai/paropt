from sympy import Symbol
from sympy.matrices import *
from sympy.printing import *
from sympy import sin, cos, Function, diff
from sympy.parsing import Maxima

x=Symbol('x')
y=Symbol('y')
f = Function('f')
M=Matrix([[x**2,y,0],[0,0,0],[2,2,2]])

DM=diff(M,x)

print ()
print(DM)
pprint (DM)

m = Maxima()
m.run_command("factor(8);")
m.run_command("factor(x^2 + 2*x*y + y^2);")

quit()


"""
import socket

HOST = ''                 # Symbolic name meaning the local host
PORT = 5007              # Arbitrary non-privileged port
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((HOST, PORT))
s.listen(1)
conn, addr = s.accept()
print ('Connected by', addr)
while 1:

    # import pudb; pu.db
    data = conn.recv(1024)
    print ("Received:", data)
    if not data: break
    conn.send(data)
conn.close()
"""
