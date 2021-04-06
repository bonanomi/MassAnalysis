# MassAnalysis
To be updated.
Set up and run in a CMSSW working area.
[Combine tool]() required.

To run the framework (for now only setup only for 1D m4l mass fit):

	python makeDCsandWSs.py -i SM_inputs_13TeV -a sm13_1D_reco_2p7fb_CB -b -d 0 -e 0 -j 0


To create workspace(s) for the fit:

	text2workspace.py <input_card>.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass -o <output_name>.root

To run the fit:

	combine -n test -M MultiDimFit <workspace_name>.root -m 125.0 -P MH --floatOtherPOIs=1 --expectSignal=1 -t -1 --X-rtd TMCSO_PseudoAsimov=1000 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --saveWorkspace --saveToys --robustFit 1 --stepSize 0.01 -s 123456 --algo=grid --points 100
