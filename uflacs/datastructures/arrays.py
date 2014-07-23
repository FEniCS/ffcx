
import numpy


def int_array(size):
    return numpy.zeros(size, dtype=int)


def object_array(size):
    return numpy.empty(size, dtype=object)


def bool_array(size):
    return numpy.zeros(size, dtype=numpy.bool8)
