
import numpy

def sufficient_int_type(maxvalue):
    if maxvalue < 2 ** 7:
        dtype = numpy.int8
    elif maxvalue < 2 ** 15:
        dtype = numpy.int16
    elif maxvalue < 2 ** 31:
        dtype = numpy.int32
    else:
        dtype = numpy.int64
    return dtype

def sufficient_uint_type(maxvalue):
    if maxvalue < 2 ** 8:
        dtype = numpy.uint8
    elif maxvalue < 2 ** 16:
        dtype = numpy.uint16
    elif maxvalue < 2 ** 32:
        dtype = numpy.uint32
    else:
        dtype = numpy.uint64
    return dtype
