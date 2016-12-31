import os
import sys

f=open("log","r")
file1={}
for line in f:
    tokenized=line.split()
    q1=tokenized[0].split("(")[1].split(",")[0]
    q2=tokenized[1].split(")")[0]
    file1[(q1,q2)]=tokenized[2]


f=open("log2","r")
file2={}
for line in f:
    tokenized=line.split()
    if len(tokenized)<7:continue
    q1=tokenized[1].split("-")[0]
    q2=tokenized[2].split("-")[0]
    file2[(str(int(q1)),str(int(q2)))]=tokenized[6]
print(len(file1),len(file2),len(file1)-len(file2))
elems=[];error=0.0
for i,val in file1.items():
    if i in file2:
        ei=float(val)-float(file2[i])
        #if abs(ei)>0.001:print(i,ei)
        error+=ei
    else:
        elems.append(i)
        print(i,val)
        continue
print(error)
print(len(elems))
