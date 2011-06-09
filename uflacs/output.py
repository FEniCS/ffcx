
from log import error

def uflacs_assert(cond, msg):
    if not cond:
        error(msg)
