import sys, os, pwd, commands
import optparse, shlex, re
import math
# from ROOT import *
# import ROOT
from array import array
from datacardClass import *
from inputReader import *

def testDCClass():
    myClass = datacardClass()
    print('Class created')
    myClass.printSomething()
    print(myClass.ID_4mu)

    myReader2e2mu = inputReader("SM_inputs_13TeV/inputs_2e2mu.txt")
    myReader2e2mu.readInputs()
    theInputs2e2mu = myReader2e2mu.getInputs()
    myClass.theInputProp(theInputs2e2mu)
    
# def testSysClass():
#     systematics = systematicsClass( self.mH, False, self.F, theInputs)
#     systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

testDCClass()