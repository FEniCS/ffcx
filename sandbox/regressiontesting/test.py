a = eval(open('a.output').read())
b = eval(open('b.output').read())
from recdiff import *
print_recdiff(recdiff(a,b))
