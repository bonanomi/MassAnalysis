#shape :  "RooChebyChev::ggZZ(mass4l,chebPol1,chebPol2,chebPol3)"
#setup : &setup
#2e2mu      :
#VBFtagged  :
#    <<       : *setup

shape = {

    'chebPol1_3_1' : 0.918273,
    'chebPol2_3_1' : 0.131573,
    'chebPol3_3_1' : 0.0546033,

#4e        :
#VBFtagged :
#    <<       : *setup
    'chebPol1_2_1' : 0.756625,
    'chebPol2_2_1' : 0.11569,
    'chebPol3_2_1' : -0.0300415,

#4mu       :
#VBFtagged :
#    <<       : *setup
    'chebPol1_1_1' : 0.355576,
    'chebPol2_1_1' : -0.019918,
    'chebPol3_1_1' : -0.0463538

}
