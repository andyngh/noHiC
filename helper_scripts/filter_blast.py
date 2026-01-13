import sys

for l in open(sys.argv[1]):
    l_arr=l.rstrip().split("\t")
    if(float(l_arr[2])>=90.0 and float(l_arr[12])>=90.0):
        print(l.rstrip())
