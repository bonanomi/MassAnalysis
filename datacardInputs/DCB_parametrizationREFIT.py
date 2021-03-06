#---

#         _\|/_
#         (o o)
# +----oOO-{_}-OOo-----------------------------------------------------------------+
# |Double Crystal-ball parameters expressed as a function of mH around mH : 125 GeV|
# +--------------------------------------------------------------------------------+
#
#setup: &setup
#    parametrization_name: 'Linear [120-130]'
#    color: 4

shape = {

#2e2mu_inclusive:
#    <<     : *setup
    'mean_3' : '124.812245508+(0.992098516531)*(@0-125)',
    'sigma_3' : '1.34318471467+(0.00446481560236)*(@0-125)',
    'alpha_3' : '1.00699495654',
    'n_3' : '2.82780783836+(-0.0224567979969)*(@0-125)',
    'alpha2_3' : '1.40516505504',
    'n2_3' : '4.3401062829+(0.0708733850487)*(@0-125)',

#4mu:
#    <<     : *setup

    'mean_1' : '124.765108592+(0.993385375758)*(@0-125)',
    'sigma_1' : '1.06432767678+(0.00656994459132)*(@0-125)',
    'alpha_1' : '1.28661471089',
    'n_1' : '1.9579637601+(-0.00601795782582)*(@0-125)',
    'alpha2_1' : '1.89341469959',
    'n2_1' : '2.7822938539+(0.0187407842516)*(@0-125)',

#4e:
#    <<     : *setup

    'mean_2' : '124.74906717+(0.988367966883)*(@0-125)',
    'sigma_2' : '1.76135801645+(0.00163621431441)*(@0-125)',
    'alpha_2' : '0.887911007092',
    'n_2' : '4.26437487431+(-0.0529211812262)*(@0-125)',
    'alpha2_2' : '1.30606463691',
    'n2_2' : '6.2445291559+(0.247379341522)*(@0-125)'


}
