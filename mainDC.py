#! /usr/bin/env python
from scipy.special import erf
import ROOT
from array import array

import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re
import math
import time
from decimal import *
import json
from ROOT import *

from systematicsClass import *
from inputReader import *

sys.path.append('./datacardInputs')

from DCB_parametrization  import *
from DCB_parametrizationREFIT  import *
from signal_shape_parametrization_13TeV_WH_VBFtagged import *
from signal_shape_parametrization_13TeV_WH_Untagged import *
from signal_shape_parametrization_13TeV_WH import *

from signal_shape_parametrization_13TeV_ZH_VBFtagged import *
from signal_shape_parametrization_13TeV_ZH_Untagged import *
from signal_shape_parametrization_13TeV_ZH import *

from bkg_shape_parametriztion_13TeV_qqZZ import *
from bkg_shape_parametriztion_13TeV_qqZZ_Untagged import *
from bkg_shape_parametriztion_13TeV_qqZZ_VBFtagged import *

from bkg_shape_parametriztion_13TeV_ggZZ import *
from bkg_shape_parametriztion_13TeV_ggZZ_Untagged import *
from bkg_shape_parametriztion_13TeV_ggZZ_VBFtagged import *

#from DCB_parametrization_refit  import *

## ------------------------------------
##  card and workspace class
## ------------------------------------

class datacardClass:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True

    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        #ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        #ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")

    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
   
 
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theis2D, theOutputDir, theInputs,theTemplateDir="templates2D", theIncludingError=False, theVBF=False, theVBFcat=0, theREFIT=False):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = theMH
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.is2D = theis2D
        self.outputDir = theOutputDir
        self.sigMorph = False#theInputs['useCMS_zz4l_sigMELA']
        self.bkgMorph = True#theInputs['useCMS_zz4l_bkgMELA']
        self.templateDir = theTemplateDir
	self.bIncludingError=theIncludingError
        self.bVBF = theVBF
        self.VBFcat = theVBFcat

        self.isREFIT = theREFIT

        self.all_chan = theInputs['all']
        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']
        self.ttbar_chan = theInputs['ttbar']
        
        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        self.low_M = 105.0
        self.high_M = 140.0
       
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"

        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)
        systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        bins = 700
        #if(self.bIncludingError) :
        #  bins = 1000

        CMS_zz4l_mass_name = "CMS_zz4l_mass"
            
        CMS_zz4l_mass = ROOT.RooRealVar(CMS_zz4l_mass_name,CMS_zz4l_mass_name,self.low_M,self.high_M)    
        CMS_zz4l_mass.setBins(bins)

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)

	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB
        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        sigma_CB_d = 0.0
    
        ## -------- Variable Definitions -------- ##
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_e_err_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_e_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_err",float(theInputs['CMS_zz4l_mean_e_sig']))
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_e_sig",0.0,-0.99,0.99)

        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_mean_m_err_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_m_err = ROOT.RooRealVar(name,"CMS_zz4l_mean_m_err",float(theInputs['CMS_zz4l_mean_m_sig']))
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",0.0,-0.99,0.99)        

        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",0.0,-0.99,0.99)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",0.0,-0.99,0.99)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",0.0,-0.99,0.99)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",0.0,-0.99,0.99)
    
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_mean_e_err.setConstant(kTRUE)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_mean_m_err.setConstant(kTRUE)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
        CMS_zz4l_alpha2.setVal(0)
        CMS_zz4l_n2.setVal(0)

        print "signal shape nuisances"
        print "mean_e_sig ", CMS_zz4l_mean_e_sig.getVal()
        print "mean_e_err ", CMS_zz4l_mean_e_err.getVal()
        print "sigma_e ", CMS_zz4l_sigma_e_sig.getVal()
        print "mean_m_sig ",CMS_zz4l_mean_m_sig.getVal()
        print "mean_m_err ", CMS_zz4l_mean_m_err.getVal()
        print "sigma_m ", CMS_zz4l_sigma_m_sig.getVal()
        print "alpha ", CMS_zz4l_alpha.getVal()
        print "n ", CMS_zz4l_n.getVal()
        print "alpha2 ", CMS_zz4l_alpha2.getVal()
        print "n2 ", CMS_zz4l_n2.getVal()

        ## -------------------- RooFormulaVar's -------------------- ##
        rfv_n_CB = ROOT.RooFormulaVar()
        rfv_alpha_CB = ROOT.RooFormulaVar()
        rfv_n2_CB = ROOT.RooFormulaVar()
        rfv_alpha2_CB = ROOT.RooFormulaVar()
        rfv_sigma_CB = ROOT.RooFormulaVar()
        rfv_mean_CB = ROOT.RooFormulaVar()

        libSigShape = 'DCB_parametrization'
        if(self.isREFIT) :
          libSigShape+='REFIT'
        print 'signal shape from',libSigShape

        _temp = __import__(libSigShape, globals(), locals(), ['shape'], -1)
        sig = _temp.shape

        ## nonresonance shapes
        _temp = __import__('signal_shape_parametrization_13TeV_WH', globals(), locals(), ['shape'], -1)
        WH = _temp.shape
        _temp = __import__('signal_shape_parametrization_13TeV_WH_Untagged', globals(), locals(), ['shape'], -1)
        WH_untagged = _temp.shape
        _temp = __import__('signal_shape_parametrization_13TeV_WH_VBFtagged', globals(), locals(), ['shape'], -1)
        WH_VBFtagged = _temp.shape

        _temp = __import__('signal_shape_parametrization_13TeV_ZH', globals(), locals(), ['shape'], -1)
        ZH = _temp.shape
        _temp = __import__('signal_shape_parametrization_13TeV_ZH_Untagged', globals(), locals(), ['shape'], -1)
        ZH_untagged = _temp.shape
        _temp = __import__('signal_shape_parametrization_13TeV_ZH_VBFtagged', globals(), locals(), ['shape'], -1)
        ZH_VBFtagged = _temp.shape

        sigWH = { }
        sigZH = { }

        if not self.bVBF:
           sigWH = WH
           sigZH = ZH
        if (self.bVBF and self.VBFcat == 0):
           sigWH = WH_untagged
           sigZH = ZH_untagged
        if (self.bVBF and self.VBFcat == 1):
           sigWH = WH_VBFtagged
           sigZH = ZH_VBFtagged

        name = ""
        if not self.bVBF:
            name = "n_{0:.0f}".format(self.channel)
            rfv_n_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))

            name = "alpha_{0:.0f}".format(self.channel)
            rfv_alpha_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")", ROOT.RooArgList(self.MH))
           
            name = "n2_{0:.0f}".format(self.channel)
            rfv_n2_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")",ROOT.RooArgList(self.MH))

            name = "alpha2_{0:.0f}".format(self.channel)
            print name
            rfv_alpha2_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")", ROOT.RooArgList(self.MH))

            name = "WHp0_{0:.0f}".format(self.channel)
            print name
            rfv_WHp0 = ROOT.RooFormulaVar(name,"("+sigWH[name]+")", ROOT.RooArgList(self.MH))

            name = "WHp1_{0:.0f}".format(self.channel)
            rfv_WHp1 = ROOT.RooFormulaVar(name,"("+sigWH[name]+")", ROOT.RooArgList(self.MH))

            name = "WHfrac_{0:.0f}".format(self.channel)
            rfv_WHfrac = ROOT.RooFormulaVar(name,"("+sigWH[name]+")", ROOT.RooArgList(self.MH))
            print 'WH res fraction ',rfv_WHfrac.getVal()

            name = "ZHp0_{0:.0f}".format(self.channel)
            rfv_ZHp0 = ROOT.RooFormulaVar(name,"("+sigZH[name]+")", ROOT.RooArgList(self.MH))

            name = "ZHp1_{0:.0f}".format(self.channel)
            rfv_ZHp1 = ROOT.RooFormulaVar(name,"("+sigZH[name]+")", ROOT.RooArgList(self.MH))

            name = "ZHfrac_{0:.0f}".format(self.channel)
            rfv_ZHfrac = ROOT.RooFormulaVar(name,"("+sigZH[name]+")", ROOT.RooArgList(self.MH))

        else:
            name = "n_{0:.0f}".format(self.channel)
            rfv_n_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
            name = "alpha_{0:.0f}".format(self.channel)
            rfv_alpha_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")", ROOT.RooArgList(self.MH))
            name = "n2_{0:.0f}".format(self.channel)
            rfv_n2_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")",ROOT.RooArgList(self.MH))
            name = "alpha2_{0:.0f}".format(self.channel)
            rfv_alpha2_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")", ROOT.RooArgList(self.MH))
            
            name = "WHp0_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            rfv_WHp0 = ROOT.RooFormulaVar(name,"("+sigWH[name]+")", ROOT.RooArgList(self.MH))
            name = "WHp1_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            rfv_WHp1 = ROOT.RooFormulaVar(name,"("+sigWH[name]+")", ROOT.RooArgList(self.MH))
            name = "WHfrac_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            rfv_WHfrac = ROOT.RooFormulaVar(name,"("+sigWH[name]+")", ROOT.RooArgList(self.MH))

            name = "ZHp0_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            rfv_ZHp0 = ROOT.RooFormulaVar(name,"("+sigZH[name]+")", ROOT.RooArgList(self.MH))
            name = "ZHp1_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            rfv_ZHp1 = ROOT.RooFormulaVar(name,"("+sigZH[name]+")", ROOT.RooArgList(self.MH))
            name = "ZHfrac_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            rfv_ZHfrac = ROOT.RooFormulaVar(name,"("+sigZH[name]+")", ROOT.RooArgList(self.MH))

        name = "mean_{0:.0f}".format(self.channel)            
        ### mean of CB
        if (self.channel == self.ID_4mu) :
             CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar(name,"("+sig[name]+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_m_err))
        elif (self.channel == self.ID_4e) :
             CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar(name,"("+sig[name]+")"+"+@0*@1*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig,CMS_zz4l_mean_e_err))
        elif (self.channel == self.ID_2e2mu) :
             CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar(name,"("+sig[name]+")"+"+(@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig,CMS_zz4l_mean_m_err,CMS_zz4l_mean_e_err))        

        # sigma of CB 
        name = "sigma_{0:.0f}".format(self.channel)
        if (self.channel == self.ID_4mu) :
            rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+sig[name]+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))

        print "fs ", self.channel
        print "n_CB ", rfv_n_CB.getVal()
        print "alpha_CB ", rfv_alpha_CB.getVal()
        print "n2_CB ", rfv_n2_CB.getVal()
        print "alpha2_CB ", rfv_alpha2_CB.getVal()
        print "sigma_CB ", rfv_sigma_CB.getVal()
        print "mean_NoConv ", CMS_zz4l_mean_sig_NoConv.getVal()

        # bIncludingError variables
        # bIncludingError variables
        merrVarName = "CMS_zz4l_massErr"

        templateName = "{0}/Dm_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        TempFile = ROOT.TFile(templateName)
        Template = ROOT.TH2F(TempFile.Get("h_mzzDm"))
        dmBins = Template.GetYaxis().GetNbins()
        dmLow = Template.GetYaxis().GetXmin()
        dmHigh = Template.GetYaxis().GetXmax()

        MassErr = ROOT.RooRealVar(merrVarName,merrVarName,dmLow,dmHigh)
        MassErr.setBins(dmBins)
        #RelErrS = ROOT.RooFormulaVar("RelErrS","@0/@1",ROOT.RooArgList(MassErr,self.MH))
        #RelErrB = ROOT.RooFormulaVar("RelErrB","@0/@1",ROOT.RooArgList(MassErr,CMS_zz4l_mass))

        rfv_MassErr = ROOT.RooFormulaVar()
        if (self.channel == self.ID_4mu) :
            rfv_MassErr = ROOT.RooFormulaVar("rfv_MassErr_1","1.0*@0*(1+@1)",ROOT.RooArgList(MassErr,CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            rfv_MassErr = ROOT.RooFormulaVar("rfv_MassErr_2","1.0*@0*(1+@1)",ROOT.RooArgList(MassErr,CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            rfv_MassErr = ROOT.RooFormulaVar("rfv_MassErr_3","1.0*@0*TMath::Sqrt((1+@1)*(1+@2))",ROOT.RooArgList(MassErr,CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))

        
        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
   
        signalCB_ggHave = ROOT.RooDoubleCB("signalCB_ggHave","signalCB_ggHave",CMS_zz4l_mass, CMS_zz4l_mean_sig_NoConv, rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB) 

        name = "signalCB_ggH_{0:.0f}".format(self.channel) 
        signalCB_ggH = ROOT.RooDoubleCB(name,name,CMS_zz4l_mass, CMS_zz4l_mean_sig_NoConv, self.getVariable(rfv_MassErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)       
        name = "signalCB_qqH_{0:.0f}".format(self.channel)
        signalCB_VBF = ROOT.RooDoubleCB(name,name,CMS_zz4l_mass, CMS_zz4l_mean_sig_NoConv, self.getVariable(rfv_MassErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
       
        name = "resWH_{0:.0f}".format(self.channel)
        resWH = ROOT.RooDoubleCB(name,name,CMS_zz4l_mass, CMS_zz4l_mean_sig_NoConv, self.getVariable(rfv_MassErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        name = "nonresWH_{0:.0f}".format(self.channel)
        nonresWH = ROOT.RooChebychev(name,name,CMS_zz4l_mass, RooArgList(rfv_WHp0,rfv_WHp1) )
        name = "signalCB_WH_{0:.0f}".format(self.channel)
        signalCB_WH = ROOT.RooAddPdf(name,name,resWH,nonresWH,rfv_WHfrac)

        name = "resZH_{0:.0f}".format(self.channel)
        resZH = ROOT.RooDoubleCB(name,name,CMS_zz4l_mass, CMS_zz4l_mean_sig_NoConv, self.getVariable(rfv_MassErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        name = "nonresZH_{0:.0f}".format(self.channel)
        nonresZH = ROOT.RooChebychev(name,name,CMS_zz4l_mass, RooArgList(rfv_ZHp0,rfv_ZHp1) )
        name = "signalCB_ZH_{0:.0f}".format(self.channel)
        signalCB_ZH = ROOT.RooAddPdf(name,name,resZH,nonresZH,rfv_ZHfrac)

        name = "signalCB_ttH_{0:.0f}".format(self.channel)
        signalCB_ttH = ROOT.RooDoubleCB(name,name,CMS_zz4l_mass, CMS_zz4l_mean_sig_NoConv, self.getVariable(rfv_MassErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)

	#----------------------- Begin  bIncludingError PDFs
        ### sig
        print "pdferr sig"
        DmtemplateSigName = "{0}/Dm_signal_{1}.root".format(self.templateDir,self.appendName)
        DmSigTempFile = ROOT.TFile(DmtemplateSigName)
        DmSigTemplate = DmSigTempFile.Get("h_Dm")
        TemplateName = "TempDataHist_SigDm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        DmSigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(MassErr),DmSigTemplate)
        TemplateName = "HistPdf_SigDm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        name = "pdfErrS_{0:.0f}".format(self.channel)
        pdfErrS = ROOT.RooHistPdf(name,name,ROOT.RooArgSet(MassErr),DmSigTempDataHist)

        #### bkg
        ### qqzz
        print "pdferr qqZZ"
        DmtemplateBkgName = "{0}/Dm_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        DmqqzzTempFile = ROOT.TFile(DmtemplateBkgName)
        DmqqzzTemplate = DmqqzzTempFile.Get("h_mzzDm")
        TemplateName = "TempDataHist_qqzzDm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        DmqqzzTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,MassErr),DmqqzzTemplate)
        TemplateName = "HistPdf_qqzzDm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        name = "pdfErr_qqzz_{0:.0f}".format(self.channel)
        pdfErrqqZZ = ROOT.RooHistPdf(name,name,ROOT.RooArgSet(CMS_zz4l_mass,MassErr),DmqqzzTempDataHist)

        ### ggzz
        print "pdferr ggZZ"
        DmtemplateBkgName = "{0}/Dm_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        DmggzzTempFile = ROOT.TFile(DmtemplateBkgName)
        DmggzzTemplate = DmggzzTempFile.Get("h_mzzDm")
        TemplateName = "TempDataHist_ggzzDm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        DmggzzTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,MassErr),DmggzzTemplate)
        TemplateName = "HistPdf_ggzzDm_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        name = "pdfErr_ggzz_{0:.0f}".format(self.channel)
        pdfErrggZZ = ROOT.RooHistPdf(name,name,ROOT.RooArgSet(CMS_zz4l_mass,MassErr),DmggzzTempDataHist)

        ## Z+X
        print "pdferr ZX"
        name = "massErrZX_ld_sigma_{0:.0f}".format(self.channel)
        rfv_EBE_zx_ld_sigma = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_ld_sigma']+")*@0", ROOT.RooArgList(CMS_zz4l_mass))
        name = "massErrZX_ld_mean_{0:.0f}".format(self.channel)
        rfv_EBE_zx_ld_mean = ROOT.RooFormulaVar(name, "("+theInputs['relerr_zx_ld_mean']+")*@0", ROOT.RooArgList(CMS_zz4l_mass))
        name = "pdfErr_zjets_{0:.0f}".format(self.channel)
        pdfErrZX = ROOT.RooLandau(name,name, MassErr, rfv_EBE_zx_ld_mean, rfv_EBE_zx_ld_sigma)

        ##### signal EBE
        name = "sig_ggHErr_{0:.0f}".format(self.channel)
        sig_ggHErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_ggH), ROOT.RooArgSet(CMS_zz4l_mass)))
        name = "sig_qqHErr_{0:.0f}".format(self.channel)
        sig_VBFErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_VBF), ROOT.RooArgSet(CMS_zz4l_mass)))
        name = "sig_WHErr_{0:.0f}".format(self.channel)
        sig_WHErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_WH), ROOT.RooArgSet(CMS_zz4l_mass)))
        name = "sig_ZHErr_{0:.0f}".format(self.channel)
        sig_ZHErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_ZH), ROOT.RooArgSet(CMS_zz4l_mass)))
        name = "sig_ttHErr_{0:.0f}".format(self.channel)
        sig_ttHErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_ttH), ROOT.RooArgSet(CMS_zz4l_mass)))

	#------------------------------------------------end bIncludingError
            
        ## --------------------------- MELA 2D PDFS ------------------------- ##
        
        discVarName = "melaLD"
    
        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        sigTempFile = ROOT.TFile(templateSigName)
        sigTemplate = sigTempFile.Get("h_mzzD")
        sigTemplate_Up = sigTempFile.Get("h_mzzD_up")
        sigTemplate_Down = sigTempFile.Get("h_mzzD_dn")

        dBins = sigTemplate.GetYaxis().GetNbins()
        dLow = sigTemplate.GetYaxis().GetXmin()
        dHigh = sigTemplate.GetYaxis().GetXmax()
        D = ROOT.RooRealVar(discVarName,discVarName,dLow,dHigh)
        D.setBins(dBins)
        print "discVarName ", discVarName, " dLow ", dLow, " dHigh ", dHigh, " dBins ", dBins
        
        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate)
        TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Up)
        TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Down)

        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ggH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ggH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_VBF = ROOT.RooHistPdf(TemplateName,TemplateName,RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_VBF_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_VBF_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_VBF_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_VBF_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_WH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_WH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_WH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_WH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ZH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ZH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ttH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ttH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ttH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplatePdf_ttH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)

        funcList_ggH = ROOT.RooArgList()  
        funcList_VBF = ROOT.RooArgList()
        funcList_WH  = ROOT.RooArgList()
        funcList_ZH  = ROOT.RooArgList()
        funcList_ttH = ROOT.RooArgList()

        funcList_ggH.add(sigTemplatePdf_ggH)
        funcList_VBF.add(sigTemplatePdf_VBF)
        funcList_WH.add(sigTemplatePdf_WH)
        funcList_ZH.add(sigTemplatePdf_ZH)
        funcList_ttH.add(sigTemplatePdf_ttH)

        if(self.sigMorph): 
           funcList_ggH.add(sigTemplatePdf_ggH_Up)
           funcList_ggH.add(sigTemplatePdf_ggH_Down)  
           funcList_VBF.add(sigTemplatePdf_VBF_Up)
           funcList_VBF.add(sigTemplatePdf_VBF_Down)  
           funcList_WH.add(sigTemplatePdf_WH_Up)
           funcList_WH.add(sigTemplatePdf_WH_Down)  
           funcList_ZH.add(sigTemplatePdf_ZH_Up)
           funcList_ZH.add(sigTemplatePdf_ZH_Down)  
           funcList_ttH.add(sigTemplatePdf_ttH_Up)
           funcList_ttH.add(sigTemplatePdf_ttH_Down)  
                
        morphSigVarName = "CMS_zz4l_sigMELA_{0:.0f}".format(self.channel)
        alphaMorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if(self.sigMorph): alphaMorphSig.setConstant(False)
        else: alphaMorphSig.setConstant(True)

        morphVarListSig = ROOT.RooArgList()
        if(self.sigMorph): 
          morphVarListSig.add(alphaMorphSig)  ## just one morphing for all signal processes        

        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ggH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_VBF = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_VBF,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_WH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ZH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        sigTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ttH,morphVarListSig,1.0,1)

        ##### 2D -> mass4l + KD  or mass4l,delta_m + KD
        name = "sigCB2d_ggH_{0:.0f}".format(self.channel)
        sigCB2d_ggH = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(sig_ggHErr,signalCB_ggH, self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(D) ) )
        name = "sigCB2d_qqH_{0:.0f}".format(self.channel)
        sigCB2d_VBF = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(sig_VBFErr,signalCB_VBF, self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(D) ) )
        name = "sigCB2d_WH_{0:.0f}".format(self.channel)
        sigCB2d_WH = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(sig_WHErr,signalCB_WH, self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(D) ) )
        name = "sigCB2d_ZH_{0:.0f}".format(self.channel)
        sigCB2d_ZH = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(sig_ZHErr,signalCB_ZH, self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(D) ) )
        name = "sigCB2d_ttH_{0:.0f}".format(self.channel)
        sigCB2d_ttH = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(sig_ttHErr,signalCB_ttH, self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(D) ) )        

        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution

        lib_qqZZShape = 'bkg_shape_parametriztion_13TeV_qqZZ'

        _qqZZtemp = __import__(lib_qqZZShape, globals(), locals(), ['shape'], -1)
        qqZZtmp = _qqZZtemp.shape

        _qqZZ_untaggedtemp = __import__(lib_qqZZShape+'_Untagged', globals(), locals(), ['shape'], -1)
        qqZZ_untaggedtmp = _qqZZ_untaggedtemp.shape
        _qqZZ_VBFtaggedtemp = __import__(lib_qqZZShape+'_VBFtagged', globals(), locals(), ['shape'], -1)
        qqZZ_VBFtaggedtmp = _qqZZ_VBFtaggedtemp.shape
          
        ## ggZZ contribution

        lib_ggZZShape = 'bkg_shape_parametriztion_13TeV_ggZZ'
        _ggZZtemp = __import__(lib_ggZZShape, globals(), locals(), ['shape'], -1)
        ggZZtmp = _ggZZtemp.shape
        _ggZZ_untaggedtemp = __import__(lib_ggZZShape+'_Untagged', globals(), locals(), ['shape'], -1)
        ggZZ_untaggedtmp = _ggZZ_untaggedtemp.shape
        _ggZZ_VBFtaggedtemp = __import__(lib_ggZZShape+'_VBFtagged', globals(), locals(), ['shape'], -1)
        ggZZ_VBFtaggedtmp = _ggZZ_VBFtaggedtemp.shape

        #################

        qqZZ = { }
        ggZZ = { }

        if not self.bVBF:
           qqZZ = qqZZtmp
           ggZZ = ggZZtmp
        if (self.VBFcat == 0 and self.bVBF):
           qqZZ = qqZZ_untaggedtmp
           ggZZ = ggZZ_untaggedtmp
        if (self.VBFcat == 1 and self.bVBF):
           qqZZ = qqZZ_VBFtaggedtmp
           ggZZ = ggZZ_VBFtaggedtmp

        ## qqZZ contribution
        if not self.bVBF:
            name = "chebPol1_{0:.0f}".format(self.channel)
            qqzz_c0 = ROOT.RooRealVar("qqzz_"+name,"qqzz_"+name,qqZZ[name])
            name = "chebPol2_{0:.0f}".format(self.channel)
            qqzz_c1 = ROOT.RooRealVar("qqzz_"+name,"qqzz_"+name,qqZZ[name])
            name = "chebPol3_{0:.0f}".format(self.channel)
            qqzz_c2 = ROOT.RooRealVar("qqzz_"+name,"qqzz_"+name,qqZZ[name])
        else:
            name = "chebPol1_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            qqzz_c0 = ROOT.RooRealVar("qqzz_"+name,"qqzz_"+name,qqZZ[name])
            name = "chebPol2_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            qqzz_c1 = ROOT.RooRealVar("qqzz_"+name,"qqzz_"+name,qqZZ[name])
            name = "chebPol3_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            qqzz_c2 = ROOT.RooRealVar("qqzz_"+name,"qqzz_"+name,qqZZ[name])

        name = "bkg_qqzzTmp_{0:.0f}".format(self.channel)
        bkg_qqzz = ROOT.RooChebychev(name,name,CMS_zz4l_mass,RooArgList(qqzz_c0,qqzz_c1,qqzz_c2))

        ## ggZZ contribution
        if not self.bVBF:
            name = "chebPol1_{0:.0f}".format(self.channel)
            ggzz_c0 = ROOT.RooRealVar("ggzz_"+name,"ggzz_"+name,ggZZ[name])
            name = "chebPol2_{0:.0f}".format(self.channel)
            ggzz_c1 = ROOT.RooRealVar("ggzz_"+name,"ggzz_"+name,ggZZ[name])
            name = "chebPol3_{0:.0f}".format(self.channel)
            ggzz_c2 = ROOT.RooRealVar("ggzz_"+name,"ggzz_"+name,ggZZ[name])
        else:
            name = "chebPol1_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            ggzz_c0 = ROOT.RooRealVar("ggzz_"+name,"ggzz_"+name,ggZZ[name])
            name = "chebPol2_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            ggzz_c1 = ROOT.RooRealVar("ggzz_"+name,"ggzz_"+name,ggZZ[name])
            name = "chebPol3_{0:.0f}_{1}".format(self.channel,self.VBFcat)
            ggzz_c2 = ROOT.RooRealVar("ggzz_"+name,"ggzz_"+name,ggZZ[name])

        name = "bkg_ggzzTmp_{0:.0f}".format(self.channel)
        bkg_ggzz = ROOT.RooChebychev(name,name,CMS_zz4l_mass,RooArgList(ggzz_c0,ggzz_c1,ggzz_c2))

        ## Reducible backgrounds
        bkg_zjets = ROOT.RooGenericPdf();
        print "zjets mass shape"
        if not self.isREFIT:
            if (self.channel == self.ID_4mu):
                p0_zjets_4mu = ROOT.RooRealVar("p0_zjets_4mu","p0_zjets_4mu",134.6)
                p1_zjets_4mu = ROOT.RooRealVar("p1_zjets_4mu","p1_zjets_4mu",24.4)
                landau_zjets_4mu = ROOT.RooFormulaVar("landau_zjets_4mu","TMath::Landau(@0,@1,@2)",RooArgList(CMS_zz4l_mass,p0_zjets_4mu,p1_zjets_4mu))
                bkg_zjets = ROOT.RooGenericPdf("bkg_zjetsTMP4mu","landau_zjets_4mu", ROOT.RooArgList(landau_zjets_4mu) )
                #TF1 *f4muC = new TF1("f4muC","landau( 0 )",70,800);
                #f4muC->SetParameters(0.04276, 134.6, 24.4);

            if (self.channel == self.ID_4e):
                p0_zjets_4e = ROOT.RooRealVar("p0_zjets_4e","p0_zjets_4e",151.2)
                p1_zjets_4e = ROOT.RooRealVar("p1_zjets_4e","p1_zjets_4e",36.6)
                p2_zjets_4e = ROOT.RooRealVar("p2_zjets_4e","p2_zjets_4e",7.06)
                p3_zjets_4e = ROOT.RooRealVar("p3_zjets_4e","p3_zjets_4e",-0.00497)
                landau_zjets_4e = ROOT.RooFormulaVar("landau_zjets_4e","TMath::Landau(@0,@1,@2)",ROOT.RooArgList(CMS_zz4l_mass,p0_zjets_4e,p1_zjets_4e))
                bkg_zjets = ROOT.RooGenericPdf("bkg_zjetsTMP4e","landau_zjets_4e*(1+TMath::Exp(p2_zjets_4e+p3_zjets_4e*CMS_zz4l_mass))", ROOT.RooArgList(landau_zjets_4e, p2_zjets_4e, p3_zjets_4e,CMS_zz4l_mass) )
                #TF1 *f4eC = new TF1("f4eC", "landau( 0 )*(1 + exp( pol1(3))) + [5]*(TMath::Landau(x, [6], [7]))", 70, 800);
                #f4eC->SetParameters(4.404e-05, 151.2, 36.6, 7.06, -0.00497, 0.01446, 157.3, 26);

            if (self.channel == self.ID_2e2mu):
                p0_zjets_2e2mu = ROOT.RooRealVar("p0_zjets_2e2mu","p0_zjets_2e2mu",144.5)
                p1_zjets_2e2mu = ROOT.RooRealVar("p1_zjets_2e2mu","p1_zjets_2e2mu",25.3)
                landau_zjets_2e2mu = ROOT.RooFormulaVar("landau_zjets_2e2mu","TMath::Landau(@0,@1,@2)",RooArgList(CMS_zz4l_mass,p0_zjets_2e2mu,p1_zjets_2e2mu))
                bkg_zjets = ROOT.RooGenericPdf("bkg_zjetsTMP2e2mu","landau_zjets_2e2mu", ROOT.RooArgList(landau_zjets_2e2mu))
                #TF1 *f2e2muC = new TF1("f2e2muC","landau(0)",70,800);
                #f2e2muC->SetParameters(0.0413, 144.5, 25.3);

        if self.isREFIT:
            if (self.channel == self.ID_4mu):
                p0_zjets_4mu = ROOT.RooRealVar("p0_zjets_4mu","p0_zjets_4mu",134.1)
                p1_zjets_4mu = ROOT.RooRealVar("p1_zjets_4mu","p1_zjets_4mu",21.01)
                landau_zjets_4mu = ROOT.RooFormulaVar("landau_zjets_4mu","TMath::Landau(@0,@1,@2)",RooArgList(CMS_zz4l_mass,p0_zjets_4mu,p1_zjets_4mu))                
                bkg_zjets = ROOT.RooGenericPdf("bkg_zjetsTMP4mu","landau_zjets_4mu", ROOT.RooArgList(landau_zjets_4mu) )
                #// shapes - 4mu
                #TF1 * fm4l_4mu = new TF1("fm4l_4mu", "[0]*(TMath::Landau(x, [1], [2]))", 70, 800);
                #fm4l_4mu->SetParameters(0.04931, 134.1, 21.01);

            if (self.channel == self.ID_4e):
                p0_zjets_4e = ROOT.RooRealVar("p0_zjets_4e","p0_zjets_4e",137)
                p1_zjets_4e = ROOT.RooRealVar("p1_zjets_4e","p1_zjets_4e",27.36)
                p2_zjets_4e = ROOT.RooRealVar("p2_zjets_4e","p2_zjets_4e",13.26)
                p3_zjets_4e = ROOT.RooRealVar("p3_zjets_4e","p3_zjets_4e",-0.002485)
                landau_zjets_4e = ROOT.RooFormulaVar("landau_zjets_4e","TMath::Landau(@0,@1,@2)",ROOT.RooArgList(CMS_zz4l_mass,p0_zjets_4e,p1_zjets_4e))
                bkg_zjets = ROOT.RooGenericPdf("bkg_zjetsTMP4e","landau_zjets_4e*(1+TMath::Exp(p2_zjets_4e+p3_zjets_4e*@3))", ROOT.RooArgList(landau_zjets_4e, p2_zjets_4e, p3_zjets_4e, CMS_zz4l_mass) )
               #// shapes - 4e
               #TF1 *fm4l_4e = new TF1("fm4l_4e", "landau(0)*(1+exp(pol1(3)))",50,800);
               #fm4l_4e->SetParameters(1.083e-07, 137, 27.36, 13.26, -0.002485);

            if (self.channel == self.ID_2e2mu):        
                p0_zjets_2e2mu = ROOT.RooRealVar("p0_zjets_2e2mu","p0_zjets_2e2mu",142.8)
                p1_zjets_2e2mu = ROOT.RooRealVar("p1_zjets_2e2mu","p1_zjets_2e2mu",23.57)
                landau_zjets_2e2mu = ROOT.RooFormulaVar("landau_zjets_2e2mu","TMath::Landau(@0,@1,@2)",RooArgList(CMS_zz4l_mass,p0_zjets_2e2mu,p1_zjets_2e2mu))
                bkg_zjets = ROOT.RooGenericPdf("bkg_zjetsTMP2e2mu","landau_zjets_2e2mu", ROOT.RooArgList(landau_zjets_2e2mu))
               #// shapes - 2e2mu
               #TF1 * fm4l_2e2mu = new TF1("fm4l_2e2mu", "landau( 0 )",70,800);
               #fm4l_2e2mu->SetParameters(0.04418, 142.8, 23.57);

        ### per-event error  
        name = "bkg_qqzzErr_{0:.0f}".format(self.channel)
	bkg_qqzzErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(bkg_qqzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrqqZZ), ROOT.RooArgSet(MassErr)));
        name = "bkg_ggzzErr_{0:.0f}".format(self.channel)
	bkg_ggzzErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(bkg_ggzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrggZZ), ROOT.RooArgSet(MassErr)));
        name = "bkg_zjetsErr_{0:.0f}".format(self.channel)
	bkg_zjetsErr = ROOT.RooProdPdf(name,name, ROOT.RooArgSet(bkg_zjets), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrZX), ROOT.RooArgSet(MassErr)));


      ## ----------------- 2D BACKGROUND SHAPES --------------- ##
        templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        bkgTempFile = ROOT.TFile(templateBkgName)
        bkgTemplate = bkgTempFile.Get("h_mzzD")
        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate)
        
        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        ggbkgTempFile = ROOT.TFile(templateggBkgName)
        ggbkgTemplate = ggbkgTempFile.Get("h_mzzD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),ggbkgTemplate)
      
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist)
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),ggbkgTempDataHist)

        ##### ZX KD
        templatezxBkgName = "{0}/Dbackground_ZX_{1}.root".format(self.templateDir ,self.appendName)
        print templatezxBkgName, "file used for ZX"

        zxbkgTempFile = ROOT.TFile(templatezxBkgName)
        zxbkgTemplate = zxbkgTempFile.Get("h_mzzD")
        zxbkgTemplate_Up = zxbkgTempFile.Get("h_mzzD_up")
        zxbkgTemplate_Down = zxbkgTempFile.Get("h_mzzD_dn")

        TemplateName = "zjetsTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        zxbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),zxbkgTemplate)
        TemplateName = "zxbkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        zxbkgTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),zxbkgTemplate_Up)
        TemplateName = "zxbkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        zxbkgTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),zxbkgTemplate_Down)
        
        TemplateName = "zxbkgTemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),zxbkgTempDataHist)
        TemplateName = "zxbkgTemplatePdf_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),zxbkgTempDataHist_Up)
        TemplateName = "zxbkgTemplatePdf_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),zxbkgTempDataHist_Down)
        
        funcList_zjets = ROOT.RooArgList()  
        morphBkgVarName = "CMS_zz4l_bkgMELA"    
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()
        
        if(self.bkgMorph):
            funcList_zjets.add(bkgTemplatePdf_zjets)
            funcList_zjets.add(bkgTemplatePdf_zjets_Up)
            funcList_zjets.add(bkgTemplatePdf_zjets_Down)  
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)  
        else:
            funcList_zjets.add(bkgTemplatePdf_zjets)
            alphaMorphBkg.setConstant(True)

        TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)

        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)

        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.bVBF):
            TemplateName += "_{0}".format(self.VBFcat)
        bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_zjets,morphVarListBkg,1.0,1)

        #### bkg 2D : mass4l + KD; mass4l, delta_m + KD
        name = "bkg2d_qqzz"
        bkg2d_qqzz = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(bkg_qqzzErr,bkg_qqzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(D) ) )
        name = "bkg2d_ggzz"
        bkg2d_ggzz = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(bkg_ggzzErr,bkg_ggzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(D) ) )
        name = "bkg2d_zjets"
        bkg2d_zjets = ROOT.RooProdPdf(name,name,ROOT.RooArgSet(self.getVariable(bkg_zjetsErr,bkg_zjets,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(D) ) )

        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##
        if not self.bVBF:
            canv_name = "czz_{0}_{1}".format(self.mH,self.appendName)
        else:
            canv_name = "czz_{0}_{1}_{2}".format(self.mH,self.appendName,self.VBFcat)
        czz = ROOT.TCanvas( canv_name, canv_name, 750, 700 )
        czz.cd()
        zzframe_s = CMS_zz4l_mass.frame(45)

        ###########
        rdsprod = ROOT.RooDataSet()
        rdsprod = sig_ggHErr.generate(ROOT.RooArgSet(CMS_zz4l_mass,MassErr), 10000)
        rdsprod.plotOn(zzframe_s)
        sig_ggHErr.plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(8) )

        # signalCB_ggH
        signalCB_ggHave.plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        nonresWH.plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(2) )
        nonresZH.plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(3) )
        bkg_qqzz.plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(4) )
        bkg_ggzz.plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(5) )
        bkg_zjets.plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )

        zzframe_s.Draw()
        if not self.bVBF:
            figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
        else:
            figName = "{0}/figs/mzz_{1}_{2}_{3}.png".format(self.outputDir, self.mH, self.appendName,self.VBFcat)
        czz.SaveAs(figName)
        del czz

        if not self.bVBF:
            canv_name = "czzerr_{0}_{1}".format(self.mH,self.appendName)
        else:
            canv_name = "czzerr_{0}_{1}_{2}".format(self.mH,self.appendName,self.VBFcat)
        czzerr = ROOT.TCanvas( canv_name, canv_name, 750, 700 )
        czzerr.cd()
        zzerrframe_s = MassErr.frame(50)

        pdfErrS.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        pdfErrqqZZ.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(2) )
        pdfErrggZZ.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(3) )
        pdfErrZX.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )

        sig_ggHErr.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(1) )
        bkg_qqzzErr.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(2) )
        bkg_ggzzErr.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(3) )
        bkg_zjetsErr.plotOn(zzerrframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(4) )

        zzerrframe_s.Draw()
        if not self.bVBF:
            figName = "{0}/figs/mzzerr_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
        else:
            figName = "{0}/figs/mzzerr_{1}_{2}_{3}.png".format(self.outputDir, self.mH, self.appendName,self.VBFcat)
        czzerr.SaveAs(figName)
        del czzerr        

        ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        #CMS_zz4l_mass.setRange("shape",self.low_M,self.high_M)
        #CMS_zz4l_mass.setRange("fullrangesignal",self.low_M,self.high_M)
        #CMS_zz4l_mass.setRange("fullrange",100,1400)

        sigRate_ggH_Shape = 1
        sigRate_VBF_Shape = 1
        sigRate_WH_Shape = 1
        sigRate_ZH_Shape = 1
        sigRate_ttH_Shape = 1
             
        fs=''
        if (self.channel == self.ID_4mu) :
           fs='4mu'
        if (self.channel == self.ID_4e) :
           fs='4e'
        if (self.channel == self.ID_2e2mu) :
           fs='2e2mu'

        print 'fs is ', fs

        ## qqZZ contribution

        print 'load signal yield parameterization'  

        _ggHtemp = __import__('ggH_norm', globals(), locals(), ['ggH_norm'], -1)
        ggH_norm = _ggHtemp.ggH_norm

        _VBFtemp = __import__('VBF_norm', globals(), locals(), ['VBF_norm'], -1)
        VBF_norm = _VBFtemp.VBF_norm

        _WpHtemp = __import__('WplusH_norm', globals(), locals(), ['WplusH_norm'], -1)
        WpH_norm = _WpHtemp.WplusH_norm

        _WmHtemp = __import__('WminusH_norm', globals(), locals(), ['WminusH_norm'], -1)
        WmH_norm = _WmHtemp.WminusH_norm

        _ZHtemp = __import__('ZH_norm', globals(), locals(), ['ZH_norm'], -1)
        ZH_norm = _ZHtemp.ZH_norm

        ## what go into ws

        ###########
        print 'signal rates'

        name_ggH_0=''
        name_ggH_1=''
        name_VBF_0=''
        name_VBF_1=''

        name_WH_0=''
        name_WH_1='' 
        name_ZH_0=''
        name_ZH_1=''
        name_ttH_0=''
        name_ttH_1=''

        if (self.channel == self.ID_4mu) :
            name_ggH_0='(-293.311)+(7.20283*@0)+(-0.0596786*@0*@0)+(0.000167384*@0*@0*@0)'
            name_ggH_1='(52.427)+(-1.25894*@0)+(0.0100537*@0*@0)+(-2.66812e-05*@0*@0*@0)'
            name_VBF_0='(-37.4819)+(0.913553*@0)+(-0.00746264*@0*@0)+(2.04681e-05*@0*@0*@0)'
            name_VBF_1='(16.5082)+(-0.386758*@0)+(0.00298349*@0*@0)+(-7.53735e-06*@0*@0*@0)'
        
            name_WH_0='(0.124686)+(-0.00413317*@0)+(2.74448e-05*@0*@0)'
            name_WH_1='(0.0655818)+(-0.00110395*@0)+(4.68736e-06*@0*@0)' 
            name_ZH_0='(0.16402)+(-0.00427212*@0)+(2.57103e-05*@0*@0)'
            name_ZH_1='(0.0132879)+(-0.000271454*@0)+(1.36914e-06*@0*@0)'
            name_ttH_0='(-0.148612)+(0.0013099*@0)'
            name_ttH_1='(-0.0112657)+(9.91373e-05*@0)'
        if (self.channel == self.ID_4e) :
            name_ggH_0='(191.477)+(-4.44232*@0)+(0.0338053*@0*@0)+(-8.3743e-05*@0*@0*@0)'
            name_ggH_1='(12.3704)+(-0.295503*@0)+(0.00233948*@0*@0)+(-6.12639e-06*@0*@0*@0)'
            name_VBF_0='(67.3174)+(-1.61452*@0)+(0.0128716*@0*@0)+(-3.40853e-05*@0*@0*@0)'
            name_VBF_1='(20.8218)+(-0.496542*@0)+(0.0039245*@0*@0)+(-1.02604e-05*@0*@0*@0)'
           
            name_WH_0='(0.156858)+(-0.00391655*@0)+(2.26816e-05*@0*@0)'
            name_WH_1='(-0.0248327)+(0.00036731*@0)+(-1.315e-06*@0*@0)'
            name_ZH_0='(0.392606)+(-0.00739397*@0)+(3.51648e-05*@0*@0)'
            name_ZH_1='(-0.00507873)+(5.29367e-05*@0)+(-7.26824e-08*@0*@0)'
            name_ttH_0='(-0.0936056)+(0.000822492*@0)'
            name_ttH_1='(-0.00561134)+(5.05038e-05*@0)'
        if (self.channel == self.ID_2e2mu) :
            name_ggH_0='(449.222)+(-10.5901*@0)+(0.0819684*@0*@0)+(-0.000206978*@0*@0*@0)'
            name_ggH_1='(89.5752)+(-2.14267*@0)+(0.0170459*@0*@0)+(-4.50737e-05*@0*@0*@0)'
            name_VBF_0='(-1.34813)+(0.0518009*@0)+(-0.000639942*@0*@0)+(2.55271e-06*@0*@0*@0)'
            name_VBF_1='(24.1087)+(-0.563716*@0)+(0.00433849*@0*@0)+(-1.09329e-05*@0*@0*@0)'
           
            name_WH_0='(0.576096)+(-0.0126016*@0)+(6.70847e-05*@0*@0)'
            name_WH_1='(0.00746677)+(-0.000201808*@0)+(1.21045e-06*@0*@0)'
            name_ZH_0='(0.325854)+(-0.00806448*@0)+(4.64559e-05*@0*@0)'
            name_ZH_1='(-0.0964998)+(0.00150981*@0)+(-5.83697e-06*@0*@0)'
            name_ttH_0='(-0.212328)+(0.00186509*@0)'
            name_ttH_1='(-0.0204792)+(0.000176347*@0)'

        name = "SigRate_ggH_0_{0:.0f}".format(self.channel)
        rfvSigRate_ggH_0 = ROOT.RooFormulaVar(name,name_ggH_0,ROOT.RooArgList(self.MH))
        name = "SigRate_ggH_1_{0:.0f}".format(self.channel)
        rfvSigRate_ggH_1 = ROOT.RooFormulaVar(name,name_ggH_1,ROOT.RooArgList(self.MH))

        name = "SigRate_VBF_0_{0:.0f}".format(self.channel)
        rfvSigRate_VBF_0 = ROOT.RooFormulaVar(name,name_VBF_0,ROOT.RooArgList(self.MH))
        name = "SigRate_VBF_1_{0:.0f}".format(self.channel)
        rfvSigRate_VBF_1 = ROOT.RooFormulaVar(name,name_VBF_1,ROOT.RooArgList(self.MH))

        name = "SigRate_WH_0_{0:.0f}".format(self.channel)
        rfvSigRate_WH_0 = ROOT.RooFormulaVar(name,name_WH_0,ROOT.RooArgList(self.MH))
        name = "SigRate_WH_1_{0:.0f}".format(self.channel)
        rfvSigRate_WH_1 = ROOT.RooFormulaVar(name,name_WH_1,ROOT.RooArgList(self.MH))

        name = "SigRate_ZH_0_{0:.0f}".format(self.channel)
        rfvSigRate_ZH_0 = ROOT.RooFormulaVar(name,name_ZH_0,ROOT.RooArgList(self.MH))
        name = "SigRate_ZH_1_{0:.0f}".format(self.channel)
        rfvSigRate_ZH_1 = ROOT.RooFormulaVar(name,name_ZH_1,ROOT.RooArgList(self.MH))

        name = "SigRate_ttH_0_{0:.0f}".format(self.channel)
        rfvSigRate_ttH_0 = ROOT.RooFormulaVar(name,name_ttH_0,ROOT.RooArgList(self.MH))
        name = "SigRate_ttH_1_{0:.0f}".format(self.channel)
        rfvSigRate_ttH_1 = ROOT.RooFormulaVar(name,name_ttH_1,ROOT.RooArgList(self.MH))

        print 'decide the final one export into the ws'

        rfvSigRate_ggH = ROOT.RooFormulaVar()
        if not self.bVBF :
           rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_hzz_norm","@0+@1",ROOT.RooArgList(rfvSigRate_ggH_0,rfvSigRate_ggH_1))
           rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_hzz_norm","@0+@1",ROOT.RooArgList(rfvSigRate_VBF_0,rfvSigRate_VBF_1))
           rfvSigRate_WH = ROOT.RooFormulaVar("WH_hzz_norm","@0+@1",ROOT.RooArgList(rfvSigRate_WH_0,rfvSigRate_WH_1))
           rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_hzz_norm","@0+@1",ROOT.RooArgList(rfvSigRate_ZH_0,rfvSigRate_ZH_1))
           rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_hzz_norm","@0+@1",ROOT.RooArgList(rfvSigRate_ttH_0,rfvSigRate_ttH_1))

           print 'ggH ', rfvSigRate_ggH.getVal(),'VBF ', rfvSigRate_VBF.getVal()
           print 'WH ', rfvSigRate_ggH.getVal(),'ZH ', rfvSigRate_ZH.getVal()

        # JES uncertainty
   
        JES = ROOT.RooRealVar("JES","JES",0,-5,5) 
        gghJES = ROOT.RooRealVar("gghJES","gghJES",0.07) 
        vbfJES = ROOT.RooRealVar("vbfJES","vbfJES",0.04)
        whJES = ROOT.RooRealVar("whJES","whJES",0.05)
        zhJES = ROOT.RooRealVar("zhJES","zhJES",0.06)
        tthJES = ROOT.RooRealVar("tthJES","tthJES",0.02)

        qqzzJES = ROOT.RooRealVar("qqzzJES","qqzzES",0.09)
        ggzzJES = ROOT.RooRealVar("ggzzJES","ggzzJES",0.08)
 

        if (self.bVBF and self.VBFcat == 0) :
            rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_hzz_norm","@0-@1*@2*@3",ROOT.RooArgList(rfvSigRate_ggH_0,rfvSigRate_ggH_1,JES,gghJES))
            rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_hzz_norm","@0-@1*@2*@3",ROOT.RooArgList(rfvSigRate_VBF_0,rfvSigRate_VBF_1,JES,vbfJES))
            rfvSigRate_WH = ROOT.RooFormulaVar("WH_hzz_norm","@0-@1*@2*@3",ROOT.RooArgList(rfvSigRate_WH_0,rfvSigRate_WH_1,JES,whJES))
            rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_hzz_norm","@0-@1*@2*@3",ROOT.RooArgList(rfvSigRate_ZH_0,rfvSigRate_ZH_1,JES,zhJES))
            rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_hzz_norm","@0-@1*@2*@3",ROOT.RooArgList(rfvSigRate_ttH_0,rfvSigRate_ttH_1,JES,tthJES))

            print 'Untagged'
            print 'ggH ', rfvSigRate_ggH.getVal(),'VBF ', rfvSigRate_VBF.getVal()
            print 'WH ', rfvSigRate_WH.getVal(),'ZH ', rfvSigRate_ZH.getVal()

        if (self.bVBF and self.VBFcat == 1) :  
            rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_hzz_norm","(@0)*(1+@1*@2)",ROOT.RooArgList(rfvSigRate_ggH_1,JES,gghJES))
            rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_hzz_norm","(@0)*(1+@1*@2)",ROOT.RooArgList(rfvSigRate_VBF_1,JES,vbfJES))
            rfvSigRate_WH = ROOT.RooFormulaVar("WH_hzz_norm","(@0)*(1+@1*@2)",ROOT.RooArgList(rfvSigRate_WH_1,JES,whJES))
            rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_hzz_norm","(@0)*(1+@1*@2)",ROOT.RooArgList(rfvSigRate_ZH_1,JES,zhJES))
            rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_hzz_norm","(@0)*(1+@1*@2)",ROOT.RooArgList(rfvSigRate_ttH_1,JES,tthJES))

            print 'VBFtagged'
            print 'ggH ', rfvSigRate_ggH.getVal(),'VBF ', rfvSigRate_VBF.getVal()
            print 'WH ', rfvSigRate_WH.getVal(),'ZH ', rfvSigRate_ZH.getVal()            


        #print 'ttH just copy ggH'
        #rfvSigRate_ttH = rfvSigRate_ggH.Clone()
        #rfvSigRate_ttH.SetName("ttH_hzz_norm")

        print 'Rates ggH ',rfvSigRate_ggH.getVal(),' VBF ',rfvSigRate_VBF.getVal()
        print 'Rates WH ',rfvSigRate_WH.getVal(),' ZH ',rfvSigRate_ZH.getVal()

        ## ----------------------- BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling

        lib_bkgnorm = 'bkg_norm';
        if(self.isREFIT) :
           lib_bkgnorm+='REFIT'

        _bkgNormtemp = __import__(lib_bkgnorm, globals(), locals(), ['bkg_norm'], -1)
        bkg_norm = _bkgNormtemp.bkg_norm

        key = "qqzz_{0:.0f}_0".format(self.channel)
        bkgRate_qqzz_untagged = bkg_norm[key]
        key = "ggzz_{0:.0f}_0".format(self.channel)
        bkgRate_ggzz_untagged = bkg_norm[key]
        key = "zjets_{0:.0f}_0".format(self.channel)
        bkgRate_zjets_untagged = bkg_norm[key]
        ##
        key = "qqzz_{0:.0f}_1".format(self.channel)
        bkgRate_qqzz_VBFtagged = bkg_norm[key]
        key = "ggzz_{0:.0f}_1".format(self.channel)
        bkgRate_ggzz_VBFtagged = bkg_norm[key]
        key = "zjets_{0:.0f}_1".format(self.channel)
        bkgRate_zjets_VBFtagged = bkg_norm[key]

        if not self.bVBF :
          bkgRate_qqzz_Shape = bkgRate_qqzz_untagged + bkgRate_qqzz_VBFtagged
          bkgRate_ggzz_Shape = bkgRate_ggzz_untagged + bkgRate_ggzz_VBFtagged
          bkgRate_zjets_Shape = bkgRate_zjets_untagged + bkgRate_zjets_VBFtagged
        if (self.bVBF and self.VBFcat == 0) :
          bkgRate_qqzz_Shape = bkgRate_qqzz_untagged 
          bkgRate_ggzz_Shape = bkgRate_ggzz_untagged 
          bkgRate_zjets_Shape = bkgRate_zjets_untagged 
        if (self.bVBF and self.VBFcat == 1) :
          bkgRate_qqzz_Shape = bkgRate_qqzz_VBFtagged 
          bkgRate_ggzz_Shape = bkgRate_ggzz_VBFtagged 
          bkgRate_zjets_Shape = bkgRate_zjets_VBFtagged       

        name = "tagRatio_qqzz_{0:.0f}".format(self.channel)
        tagRatio_qqzz = RooRealVar(name,name,bkgRate_qqzz_VBFtagged/bkgRate_qqzz_untagged)
        name = "tagRatio_ggzz_{0:.0f}".format(self.channel)
        tagRatio_ggzz = RooRealVar(name,name,bkgRate_ggzz_VBFtagged/bkgRate_ggzz_untagged)     

        if (self.bVBF and self.VBFcat == 0) :
            rfvSigRate_qqzz = ROOT.RooFormulaVar("qqzz_norm","1-@0*@1*@2",ROOT.RooArgList(JES,qqzzJES,tagRatio_qqzz))
            rfvSigRate_ggzz = ROOT.RooFormulaVar("ggzz_norm","1-@0*@1*@2",ROOT.RooArgList(JES,ggzzJES,tagRatio_ggzz))

        if (self.bVBF and self.VBFcat == 1) :
            rfvSigRate_qqzz = ROOT.RooFormulaVar("qqzz_norm","1+@0*@1",ROOT.RooArgList(JES,qqzzJES))
            rfvSigRate_ggzz = ROOT.RooFormulaVar("ggzz_norm","1+@0*@1",ROOT.RooArgList(JES,ggzzJES))

        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        if not self.bVBF:
            dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
        else:
            dataFileName = "{0}/hzz{1}_{2}_{3}.root".format(dataFileDir,self.appendName,self.lumi,self.VBFcat)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        
        if (self.is2D == 0):
            if(self.bIncludingError): data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass, MassErr))
            else: data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass))
		
        if (self.is2D == 1):
            if(self.bIncludingError):
                data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D,MassErr) )
            else: data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D) )

        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""
        name_ShapeWSXSBR = ""

        if not self.bVBF:
            if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)

            if (endsInP5): name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

            if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

            name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)
        else:
            if (endsInP5): name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)

            if (endsInP5): name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)
            else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)

            if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)
            else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV_{4}.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.VBFcat)

            name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV_{2}.input.root".format(self.appendName,self.sqrts,self.VBFcat)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        #w.importClassCode(RooqqZZPdf_v2.Class(),True)
        #w.importClassCode(RooggZZPdf_v2.Class(),True)
        w.importClassCode(RooDoubleCB.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)

        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
    
        if (self.is2D == 0):
	    if not self.bIncludingError:
                	signalCB_ggH.SetNameTitle("ggH_hzz","ggH_hzz")
                	signalCB_VBF.SetNameTitle("qqH_hzz","qqH_hzz")
                	signalCB_WH.SetNameTitle("WH_hzz","WH_hzz")
                	signalCB_ZH.SetNameTitle("ZH_hzz","ZH_hzz")
                	signalCB_ttH.SetNameTitle("ttH_hzz","ttH_hzz")
                
                	getattr(w,'import')(signalCB_ggH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_VBF, ROOT.RooFit.RecycleConflictNodes())
               		getattr(w,'import')(signalCB_WH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_ZH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_ttH, ROOT.RooFit.RecycleConflictNodes())
	    else:
                	sig_ggHErr.SetNameTitle("ggH_hzz","ggH_hzz")
                	sig_VBFErr.SetNameTitle("qqH_hzz","qqH_hzz")
                	sig_WHErr.SetNameTitle("WH_hzz","WH_hzz")
                	sig_ZHErr.SetNameTitle("ZH_hzz","ZH_hzz")
                	sig_ttHErr.SetNameTitle("ttH_hzz","ttH_hzz")
                
                	getattr(w,'import')(sig_ggHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_VBFErr, ROOT.RooFit.RecycleConflictNodes())
               		getattr(w,'import')(sig_WHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_ZHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_ttHErr, ROOT.RooFit.RecycleConflictNodes())
                
        if (self.is2D == 1):
        
                    sigCB2d_ggH.SetNameTitle("ggH_hzz","ggH_hzz")
                    sigCB2d_VBF.SetNameTitle("qqH_hzz","qqH_hzz")
                    sigCB2d_WH.SetNameTitle("WH_hzz","WH_hzz")
                    sigCB2d_ZH.SetNameTitle("ZH_hzz","ZH_hzz")
                    sigCB2d_ttH.SetNameTitle("ttH_hzz","ttH_hzz")
                
                    getattr(w,'import')(sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_VBF, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_WH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_ZH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigCB2d_ttH, ROOT.RooFit.RecycleConflictNodes()) 

        if (self.is2D == 0):
		if not self.bIncludingError:
			bkg_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
			bkg_ggzz.SetNameTitle("bkg_ggzz","bkg_ggzz")
			bkg_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
            		getattr(w,'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_ggzz, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
		else:
			bkg_qqzzErr.SetNameTitle("bkg_qqzz","bkg_qqzz")
			bkg_ggzzErr.SetNameTitle("bkg_ggzz","bkg_ggzz")
			bkg_zjetsErr.SetNameTitle("bkg_zjets","bkg_zjets")
            		getattr(w,'import')(bkg_qqzzErr, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_ggzzErr, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_zjetsErr, ROOT.RooFit.RecycleConflictNodes())
            
        if (self.is2D == 1):
                getattr(w,'import')(bkg2d_qqzz,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(bkg2d_ggzz,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(bkg2d_zjets,ROOT.RooFit.RecycleConflictNodes())

        
        getattr(w,'import')(rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_VBF, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_WH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ZH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ttH, ROOT.RooFit.RecycleConflictNodes())
            
        CMS_zz4l_mass.setRange(self.low_M,self.high_M)

        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)
        
        ## --------------------------- DATACARDS -------------------------- ##

        systematics.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape, bkgRate_zjets_Shape)
        systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape,bkgRate_zjets_Shape)

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan and not self.all_chan :  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = bkgRate_qqzz_Shape
        rates['ggZZ']  = bkgRate_ggzz_Shape
        rates['zjets'] = bkgRate_zjets_Shape
        rates['ttbar'] = 0
        rates['zbb']   = 0
        

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        if not self.bVBF:
            systematics.WriteSystematics(fo, theInputs)
            systematics.WriteShapeSystematics(fo,theInputs)
        else:
            systematics.WriteSystematics(fo, theInputs,self.VBFcat)
            systematics.WriteShapeSystematics(fo,theInputs)
        fo.close()

        ## forXSxBR
        if not self.bVBF:
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
        else:
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.bVBF)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV_{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.bVBF)
            
        fo = open( name_Shape, "wb" )
        if not self.bVBF:
            self.WriteDatacard(fo, theInputs,name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        else:
            self.WriteDatacard(fo,theInputs,name_ShapeWS2,rates,data_obs.numEntries(),self.is2D,False,"",True,self.VBFcat)
        systematics_forXSxBR.WriteSystematics(fo, theInputs)
        systematics_forXSxBR.WriteShapeSystematics(fo,theInputs)
        fo.close()


    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents,is2D,isAltCard=False,AltLabel="",bVBF=False,VBFcat=""):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)
        
        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        
        if not self.bVBF:
            file.write("bin a{0} \n".format(self.channel))
        else:
            file.write("bin a{0}_{1} \n".format(self.channel,self.VBFcat))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("bin ")        

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']

        channelName1D=['ggH_hzz','qqH_hzz','WH_hzz','ZH_hzz','ttH_hzz','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
        channelName2D=['ggH_hzz','qqH_hzz','WH_hzz','ZH_hzz','ttH_hzz','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

#            channelList=['ggH{0}'.format(AltLabel),'qqZZ','ggZZ','zjets','ttbar','zbb']
         
        for chan in channelList:
            if theInputs[chan]:
                if not self.bVBF:
                    file.write("a{0} ".format(self.channel))
                else:
                    file.write("a{0}_{1} ".format(self.channel,self.VBFcat))
            else:
                if chan.startswith("ggH") and theInputs["all"] :
                    file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0
        if not (self.is2D == 1):
            for chan in channelList:
                if theInputs[chan]:
                    file.write("{0} ".format(channelName1D[i]))
                i+=1
        else:
            for chan in channelList:
#                print 'checking if ',chan,' is in the list of to-do'
                if theInputs[chan]:
                    file.write("{0} ".format(channelName2D[i]))
 #                   print 'writing in card index=',i,'  chan=',chan
                    i+=1
                else:
                    if chan.startswith("ggH") and theInputs["all"] :
                        file.write("{0} ".format(channelName2D[i]))
  #                      print 'writing in card TOTAL SUM, index=',i,'  chan=',chan,'  ',channelName2D[i]
                        i+=1
        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        for chan in channelList:
            if theInputs[chan] or (chan.startswith("ggH") and theInputs["all"]):
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


        
    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggH']: counter+=1
        if inputs['qqH']: counter+=1
        if inputs['WH']:  counter+=1
        if inputs['ZH']:  counter+=1
        if inputs['ttH']: counter+=1
        if inputs['all']: counter+=1
        
        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1
        
        return counter

