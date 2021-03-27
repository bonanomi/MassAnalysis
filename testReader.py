import sys, os, pwd, commands
import optparse, shlex, re
import math
from array import array
from inputReader import *

#define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-fs', '--fState', dest='finState', type='string', default="4d",    help='Final State')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

#define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()

inputDir = 'SM_inputs_13TeV'
def testInput(inputDir, finState):
    fname = "/inputs_%s.txt" %finState
    myReader = inputReader(inputDir+fname)
    myReader.readInputs()
    theInputs = myReader.getInputs()

    print theInputs

finState = '4mu'
testInput(inputDir, finState)