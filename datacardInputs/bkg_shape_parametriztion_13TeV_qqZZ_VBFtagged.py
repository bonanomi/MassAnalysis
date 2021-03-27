#shape : "RooChebyChev::qqZZ(mass4l,chebPol1,chebPol2,chebPol3)"
#setup : &setup
shape = {

#2e2mu     :
#VBFtagged :
#    <<       : *setup

    'chebPol1_3_1' : 0.316454,
    'chebPol2_3_1' : -0.12121,
    'chebPol3_3_1' : 0.0692207,

#4e        :
#VBFtagged :
#    <<       : *setup
    'chebPol1_2_1' : 0.241958,
    'chebPol2_2_1' : -0.0568398,
    'chebPol3_2_1' : 0.0206944,

#4mu       :
#VBFtagged :
#    <<       : *setup
    'chebPol1_1_1' : 0.267476,
    'chebPol2_1_1' : 0.0140923,
    'chebPol3_1_1' : 0.110047

}
