
CXX=g++
RM=rm
PYTHON=python

TESTAPP=testall

GENDIR=generated
CPPDIR=cpp
SCRIPTSDIR=scripts
PYDIR=py

COMMONHEADERS=${CPPDIR}/common/*.h
TEST_PY=$(wildcard ${PYDIR}/test_*.py)
TEST_PY_BASE=$(subst ${PYDIR}/,,${TEST_PY})
TEST_PROG=${TEST_PY_BASE:.py=}
TEST_SRC=$(addsuffix .h, $(addprefix ${CPPDIR}/,${TEST_PROG}))
#TEST_SRC=$(wildcard ${CPPDIR}/test_*.h)
TEST_SRC_BASE=$(subst ${CPPDIR}/,,${TEST_SRC})
#TEST_PROG=${TEST_SRC_BASE:.h=}
TEST_RUN=$(addprefix run_,${TEST_PROG})
TEST_GC=$(addprefix ${GENDIR}/gc_,${TEST_SRC_BASE})
TEST_MAINS=$(addsuffix .cpp, $(addprefix ${GENDIR}/main_,${TEST_PROG}))

CXXFLAGS=-I${GTEST_DIR}/include -I${GENDIR} -I${CPPDIR}
LIBDIRS=-L${GTEST_DIR}/lib
LIBS=-lgtest -lpthread

.PHONY: show clean purge all run run_testall run_test_*

.SECONDARY:

all: run_testall

run: ${TEST_RUN}

run_testall: ${TESTAPP}
	./${TESTAPP}

run_test_%: test_%
	./$<

${TESTAPP}: ${CPPDIR}/main_testall.cpp ${CPPDIR}/main_testall.h ${TEST_SRC} ${TEST_GC} ${COMMONHEADERS}
	${CXX} -o $@ ${CXXFLAGS} ${LIBDIRS} $< ${LIBS}

test_%: ${GENDIR}/main_test_%.cpp ${GENDIR}/gc_test_%.h ${CPPDIR}/test_%.h ${COMMONHEADERS}
	${CXX} -o $@ ${CXXFLAGS} ${LIBDIRS} $< ${LIBS}

${CPPDIR}/test_%.h:
	${PYTHON} ${SCRIPTSDIR}/init_test_header.py $@

${GENDIR}/gc_test_%.h ${GENDIR}/main_test_%.cpp: ${PYDIR}/test_%.py ${PYDIR}/codegentestcase.py
	${PYTHON} $<

show:
	@echo
	@echo src
	@echo   ${TEST_SRC}
	@echo
	@echo src base
	@echo   ${TEST_SRC_BASE}
	@echo
	@echo prog
	@echo   ${TEST_PROG}
	@echo
	@echo run
	@echo   ${TEST_RUN}
	@echo
	@echo gc
	@echo   ${TEST_GC}
	@echo
	@echo mains
	@echo   ${TEST_MAINS}
	@echo

clean:
	${RM} -f ${GENDIR}/*
	${RM} -f `find -name \*.pyc`
	${RM} -f *.log

purge: clean
	${RM} -f ${TEST_PROG} ${TESTAPP}
	${RM} -f `find -name \*.py~`
