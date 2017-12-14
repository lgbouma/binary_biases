# "Nominal" models (this just means I arbitrarily decided Z_2 = 0.5)
python numerical_models.py --modelnumber 1 --ZsubTwo 0.5 --binaryfrac 0.1
python numerical_models.py --modelnumber 2 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
python numerical_models.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
python numerical_models.py --modelnumber 4 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 22.5
python numerical_models.py --modelnumber 5 --ZsubTwo 0.5 --binaryfrac 0.1 --upperradiuscutoff 22.5

# Check vs Z_2
python numerical_models.py --modelnumber 3 --ZsubTwo 0.0 --binaryfrac 0.44 --upperradiuscutoff 22.5
python numerical_models.py --modelnumber 3 --ZsubTwo 0.25 --binaryfrac 0.44 --upperradiuscutoff 22.5

# Check vs upper radius cutoff
python numerical_models.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 15.0
python numerical_models.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 17.5
python numerical_models.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 20.0
python numerical_models.py --modelnumber 3 --ZsubTwo 0.5 --binaryfrac 0.44 --upperradiuscutoff 25.0
