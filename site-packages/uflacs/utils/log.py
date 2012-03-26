
class UflacException(Exception):
    pass

def info(msg):
    print msg

def warning(msg):
    print msg

def error(msg):
    print msg
    raise UflacException(msg)
