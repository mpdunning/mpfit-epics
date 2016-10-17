TOP=../..

include $(TOP)/configure/CONFIG
LIBRARY_IOC += mpfit
INC += mpfit.h fitFuncs.h
mpfit_SRCS += mpfit.c fitFuncs.c
mpfit_LIBS += $(EPICS_BASE_IOC_LIBS)
include $(TOP)/configure/RULES
