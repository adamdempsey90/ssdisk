import sys                                                                                               
import os                                                                                                
here = os.path.realpath("."); path = sys.path.insert(0,here+"/scripts")                                  
try:                                                                                                     
    import test as T                                                                                     
except ImportError:                                                                                      
    print "\nError!!! test module can not be imported. Be sure that you're executing the test from the main directory, using make for that.\n"

description1 = """Running fargo setup without ghostsX on the CPU.\n"""                                   
description2 = """Running fargo setup with ghostsX on the GPU.\n"""                                      
GhostsFargoTest = T.GenericTest(testname = "fargo_Ghosts_TEST", 
                           flags1 = "SETUP=fargo GHOSTSX=0 PARALLEL=0 FARGO_DISPLAY=NONE GPU=0",         
                           launch1 = "./fargo3d -m",
                           description1 = description1,                                                  
                           flags2 = "SETUP=fargo GHOSTSX=1 PARALLEL=0 FARGO_DISPLAY=NONE GPU=1",                   
                           launch2 = "./fargo3d -D 0 -m",
                           description2 = description2,
                           parfile = "setups/fargo/fargo.par",                                           
                           verbose = False,
                           plot=False,
                           field = "gasdens",                                                            
                           compact = True,
                           parameters = {'dt':0.2, 'ninterm':1, 'ntot':1,                                
                                         'ny':64, 'nx':64, 'nz':1},
                           clean = True) 
GhostsFargoTest.run()  
