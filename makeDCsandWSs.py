#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.08.30
# by Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array
from datacardClass import *
from inputReader import *

#define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-i', '--input', dest='inputDir', type='string', default="",    help='inputs directory')
    parser.add_option('-d', '--is2D',   dest='is2D',       type='int',    default=1,     help='is2D (default:1)')
    parser.add_option('-a', '--append', dest='appendName', type='string', default="",    help='append name for cards dir')
    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    parser.add_option('-t', '--templateDir', type='string', dest='templateDir', default="templates2D" ,help='directory with 2D template histos')
    parser.add_option('-e', '--massError',   dest='massError',       type='int',    default=0,     help='massError (default:0)')
    parser.add_option('-j', '--jet', dest='useJET', type='int', default=0, help='useJET (default:0)')

    parser.add_option('-r', '--refit', dest='useREFIT', type='int', default=0, help='useREFIT (default:0)')
    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.is2D != 0 and opt.is2D != 1):
        print 'The input '+opt.is2D+' is unkown for is2D.  Please choose 0 or 1. Exiting...'
        sys.exit()

    if (opt.appendName == ''):
        print 'Please pass an append name for the cards directory! Exiting...'
        sys.exit()

    if (opt.inputDir == ''):
        print 'Please pass an input directory! Exiting...'
        sys.exit()

# define make directory function
def makeDirectory(subDirName):
    if (not os.path.exists(subDirName)):
        cmd = 'mkdir -p '+subDirName
        status, output = commands.getstatusoutput(cmd)
        if status !=0:
            print 'Error in creating submission dir '+subDirName+'. Exiting...'
            sys.exit()
    else:
        print 'Directory '+subDirName+' already exists. Exiting...'
        sys.exit()


#define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()



def creationLoop(directory):
    global opt, args

    startMass=[125.0]
    stepSizes=[0.1]
    endVal=[1]


    myClass = datacardClass()
    myClass.loadIncludes()

    if(opt.useREFIT) :
       opt.templateDir = 'templates2D_refit'

    # Always use -j 0, so leave only this
    # inputReader doesn't take any argument, so it should
    # work out of the box
    # Input dir set with -i option and will be SM_inputs_13TeV
    # Can be tested with python testReader.py
    myReader4e = inputReader(opt.inputDir+"/inputs_4e.txt")
    myReader4e.readInputs()
    theInputs4e = myReader4e.getInputs()
    
    myReader4mu = inputReader(opt.inputDir+"/inputs_4mu.txt")
    myReader4mu.readInputs()
    theInputs4mu = myReader4mu.getInputs()
    
    myReader2e2mu = inputReader(opt.inputDir+"/inputs_2e2mu.txt")
    myReader2e2mu.readInputs()
    theInputs2e2mu = myReader2e2mu.getInputs()
    
    a=0
    while (a < len(startMass) ):
	
	c = 0
        while (c < endVal[a] ): 
            
            mStart = startMass[a]
            step = stepSizes[a]
            mh = mStart + ( step * c ) 
            mhs = str(mh).replace('.0','')

            print mh
            # Always use -j 0 so leave only this
            # Create cards and workspace for the three final state
            # The only thing that changes is mh, i.e. the mass value used for the ws
            print 'useJET == 0'
            makeDirectory(directory+'/HCG/'+mhs)
            makeDirectory(directory+'/HCG_XSxBR/'+mhs)

            # For the moment is2D will be False always
            # For the moment w/o massError and w/o refit
            # This will run 3 times, one for each directory in subdir = ['HCG','HCG_XSxBR','figs']
            myClass.makeCardsWorkspaces(mh,opt.is2D,directory,theInputs4e,opt.templateDir, opt.massError,0,0,opt.useREFIT)
            myClass.makeCardsWorkspaces(mh,opt.is2D,directory,theInputs4mu,opt.templateDir,opt.massError,0,0,opt.useREFIT)
            myClass.makeCardsWorkspaces(mh,opt.is2D,directory,theInputs2e2mu,opt.templateDir,opt.massError,0,0,opt.useREFIT)

            c += 1
	a += 1

# the main procedure
def makeDCsandWSs():
    
    # parse the arguments and options
    global opt, args
    parseOptions()

    if (opt.appendName != ''):
        dirName = 'cards_'+opt.appendName
    
    # Destination directories for the outputs created by the code
    subdir = ['HCG','HCG_XSxBR','figs']

    for d in subdir:
        makeDirectory(dirName+'/'+d)

    creationLoop(dirName)

    sys.exit()



# run the create_RM_cfg() as main()
if __name__ == "__main__":
    makeDCsandWSs()


